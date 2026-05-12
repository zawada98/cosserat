# -*- coding: utf-8 -*-
"""
ctr_gui.py
==========
Tkinter-based GUI bridge for the two-tube CTR simulation.

ROLE IN THE PIPELINE
--------------------
This file is purely the *human-machine interface* layer. It owns:
  - a tk.Tk window with sliders, knobs and an Initialize button
  - a tk.Toplevel window hosting a live matplotlib plot of the contact
    profile (normal gap δn vs curvilinear abscissa of Pc_B along Tube_3)
  - a daemon thread running tkinter mainloop
  - a threading.Lock-guarded shared dict that mediates state with SOFA

It does NOT own any DOFs, ForceFields, or SOFA components. It does NOT
know about Cosserat strains or contact pairs. It only knows that there
are two tubes; each can be translated and rotated; there is a dt; and
that the SOFA side periodically posts a (curvilinear-abscissa, normal-
gap) array pair into a shared buffer for the live plot to render.

THREADING ARCHITECTURE (the part that took the most thought)
-----------------------------------------------------------
SOFA's main thread runs runSofa's GLFW loop and owns the GIL state during
animation steps. Tkinter cannot share that thread (this is the same root
cause that broke the live matplotlib popup in init_monitoring.py).
The GUI therefore runs in its own daemon thread.

Four cross-thread channels are used:

  GUI -> SOFA  (high-frequency, every Scale change)
    Slider/knob commands write into self._shared under self._lock.
    CTRController.snapshot() returns a dict copy each onAnimateBeginEvent.
    Lock-protected dict access is the only safe pattern here -- direct
    attr writes on Sofa.Core.Data fields from a non-main thread can race
    with the C++ side.

  SOFA -> GUI  (one-shot, init complete signal)
    InitializationMonitor calls bridge.signal_init_complete(), which
    only sets a flag in the dict. The GUI's own poll loop
    (root.after(POLL_INTERVAL_MS, _poll)) detects the flag and enables
    the control widgets.
    We do NOT mutate Tk widgets from the SOFA thread -- Tk is not
    thread-safe and that pattern crashes under load. "GUI thread polls a
    flag, GUI thread mutates its own widgets" is the accepted workaround.

  GUI -> SOFA  (one-shot, Initialize button press)
    The button handler sets root_node.animate = True directly. This is
    a single Data field write; the worst-case race is one extra step or
    one missed step before the change takes effect, which is benign.

  SOFA -> GUI  (high-frequency, live contact profile)
    LiveContactMonitor calls bridge.push_contact_profile(step, abscissae,
    gaps) every N simulation steps.  The arrays are copied into a single-
    frame "latest value" slot in the shared dict (no queue -- if Tk is
    slower than SOFA, intermediate frames are silently overwritten,
    which is exactly the behavior we want for a live monitor).  The Tk
    thread's _poll_plot loop reads the slot at PLOT_POLL_INTERVAL_MS,
    updates pre-built matplotlib artists with set_data / set_offsets,
    and requests a redraw via draw_idle().  No artist creation per
    frame, no GUI-toolkit calls from the SOFA side.

PHASE STATE MACHINE (mirrored in self._shared['phase'])
-------------------------------------------------------
  'waiting'      : scene loaded, root_node.animate = False (paused).
                   GUI shows the Initialize button enabled, all controls
                   greyed out.
  'initializing' : Initialize button has been clicked. animate = True;
                   inner tube relaxes from conform-to-outer.  Controls
                   stay greyed out. InitializationMonitor watches the
                   FramesMO velocity field and signals when settled.
  'control'      : InitializationMonitor has fired signal_init_complete.
                   Controls are enabled. dt is unchanged from the init
                   phase value; user can change it via the Spinbox + Apply
                   button at any time.  (Earlier versions auto-bumped dt
                   here and crashed the solver -- see signal_init_complete
                   for the post-mortem.)

LIMITATIONS
-----------
  - macOS: Tk requires the main thread; this bridge will fail there.
  - Linux: requires python3-tk (`apt install python3-tk`).
  - Closing the main control window does NOT shut down SOFA -- it just
    hides the GUI.  Closing the live plot window also only hides it;
    the SOFA-side push keeps updating the buffer.  To stop the
    simulation, close runSofa.
  - matplotlib is OPTIONAL: the control GUI works without it.  Only the
    live plot Toplevel is suppressed if matplotlib import fails.
"""

import threading
import tkinter as tk
from tkinter import ttk

import numpy as np      # [ADDED] for live-plot data buffers

# Matplotlib is optional -- the control GUI works without it; only the live
# contact-profile window depends on it.  The import is wrapped so a missing
# matplotlib does not prevent the rest of the GUI from coming up.
try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
    _MPL_OK = True
except ImportError as _e:
    print(f"[CTRGui] matplotlib not available, live plot disabled: {_e!r}")
    _MPL_OK = False


# Sentinel for shared['dt_request']: None means "no pending change".
_NO_DT_REQUEST = None


class CTRGuiBridge:
    """
    Owns the Tkinter GUI thread and the shared-state dict.

    Construct in createScene() and stash on the CTRController so Python's
    GC can't collect it while SOFA runs.
    """

    POLL_INTERVAL_MS      = 100   # how often the GUI polls shared['init_complete']
    PLOT_POLL_INTERVAL_MS = 50    # ~20 Hz repaint cadence for the live plot

    # ==================================================================
    #  Construction
    # ==================================================================
    def __init__(self,
                 root_node,
                 t1_max_translation_m=0.04,        # outer slider max [m]
                 t2_max_translation_m=0.08,        # inner slider max [m]
                 max_rotation_rate_deg_per_step=2.0,   # knob max magnitude [deg/step]
                 init_dt=1e-4,                     # dt held during init phase
                 default_control_dt=1e-3,          # dt auto-applied at control entry
                 dt_min=1e-6, dt_max=1e-1,         # allowed dt range in GUI
                 default_trans_step_um=50.0,       # default slider chase speed [µm/step]
                 trans_step_min_um=1.0,            # Spinbox lower bound  [µm/step]
                 trans_step_max_um=500.0):         # Spinbox upper bound  [µm/step]

        self._root_node               = root_node
        self._t1_max                  = float(t1_max_translation_m)
        self._t2_max                  = float(t2_max_translation_m)
        self._max_rot_rate            = float(max_rotation_rate_deg_per_step)
        self._init_dt                 = float(init_dt)
        self._default_control_dt      = float(default_control_dt)
        self._dt_min                  = float(dt_min)
        self._dt_max                  = float(dt_max)
        self._default_trans_step_um   = float(default_trans_step_um)
        self._trans_step_min_um       = float(trans_step_min_um)
        self._trans_step_max_um       = float(trans_step_max_um)

        # ---- Shared state (guarded by self._lock; readable from any thread) ----
        self._lock = threading.Lock()
        self._shared = {
            'phase':                    'waiting',
            'init_pressed':             False,
            'init_complete':            False,
            't1_translation_target_m':  0.0,
            't2_translation_target_m':  0.0,
            # [MODIFIED] Rotation knob value is now in deg PER SIMULATION STEP,
            # not deg/sec sim time.  Same logic as translation: per-step is
            # the unit that gives a dt-independent real-time perception, which
            # is what the user feels.  The previous deg/sec sim time gave
            # ~3.6 deg/real-sec at full slider for this scene's dt = 1e-4 --
            # easy to mistake for "rotation does not work".
            't1_rotation_rate_deg_per_step':   0.0,
            't2_rotation_rate_deg_per_step':   0.0,
            'dt_request':               _NO_DT_REQUEST,
            # Live-tunable translation chase speed in *meters per simulation
            # step* (NOT m/s).
            #
            # Why per-step rather than per-second-sim-time:
            #   The user's perception of motion speed is in real wall-clock
            #   time.  Wall-clock_speed = step_size * sim_steps_per_real_sec.
            #   sim_steps_per_real_sec depends on solver throughput, which
            #   is roughly dt-independent for a non-adaptive solver.  So
            #   step_size in m/step gives a real-time speed that does NOT
            #   change when the user changes dt.
            #   Per-second-sim-time, by contrast, multiplies by dt to get
            #   step_size, so doubling dt doubles the per-step jump that
            #   the constraint solver feels -- exactly the wrong coupling.
            #
            # Also: per-step is what the constraint solver "feels" each
            # step.  Capping per-step keeps kinematic shocks bounded
            # regardless of dt.
            'translation_step_m':       float(default_trans_step_um) * 1e-6,

            # ---- [ADDED] Live contact-profile buffer ----------------------
            # Single-frame "latest value" slot (no queue).  SOFA writes via
            # push_contact_profile(); GUI reads via _poll_plot() and clears
            # the dirty flag.  If the GUI repaint is slower than the SOFA
            # push rate, intermediate frames are silently overwritten --
            # exactly what we want for a live monitor: always show the
            # freshest data, never a backlog.
            'contact_profile': {
                'step':      0,
                'abscissae': None,    # numpy array, meters (or None)
                'gaps_m':    None,    # numpy array, meters (or None)
                'dirty':     False,   # set True on push, False on consume
            },
        }

        # ---- GUI-thread-only state (only touched inside _run / _on_*) ----
        self._tk_root                 = None
        self._init_button             = None
        self._status_label            = None
        self._control_frames          = []
        self._control_widgets_enabled = False
        # Tk variables
        self._var_t1_tx = self._var_t1_rot = None
        self._var_t2_tx = self._var_t2_rot = None
        self._var_dt    = None
        self._var_trans_step = None     # [MODIFIED] translation step Spinbox (µm/step)
        # Live value labels
        self._lbl_t1_tx = self._lbl_t1_rot = None
        self._lbl_t2_tx = self._lbl_t2_rot = None

        # [ADDED] Live-plot Tk Toplevel + matplotlib state (GUI-thread only)
        self._plot_window  = None
        self._plot_canvas  = None
        self._plot_ax      = None
        self._plot_line    = None       # connecting line artist (Line2D)
        self._plot_scatter = None       # per-pair markers (PathCollection)

        # Start the GUI thread. Daemon -> dies with the SOFA process.
        self._thread = threading.Thread(target=self._run,
                                        name='CTRGuiThread',
                                        daemon=True)
        self._thread.start()

    # ==================================================================
    #  SOFA-side API  (thread-safe; callable from any thread)
    # ==================================================================
    def snapshot(self):
        """Atomic dict copy of the shared state. Called every step."""
        with self._lock:
            return dict(self._shared)

    def consume_dt_request(self):
        """Atomically read AND clear the pending dt_request, or None."""
        with self._lock:
            v = self._shared['dt_request']
            self._shared['dt_request'] = _NO_DT_REQUEST
            return v

    def signal_init_complete(self):
        """
        Called by InitializationMonitor when relaxation has settled.
        Side effects (all in shared dict):
          - phase             -> 'control'
          - init_complete     -> True

        The GUI poll loop notices init_complete and enables the control
        widgets on its own thread.

        Note on dt -- and why this method does NOT change it
        ----------------------------------------------------
        Earlier versions of this method auto-queued
            dt_request = default_control_dt
        here, so the controller would replace the (slow) init dt with the
        (faster) default control dt on the very next AnimateBegin.  That
        caused the solver to explode at the init->control transition:
        at the moment this fires, the system is in constraint equilibrium
        FOR dt = init_dt, and a sudden 10x dt jump (1e-4 -> 1e-3) overruns
        the BlockGaussSeidel constraint solve (tolerance=1e-5, maxIter=500
        in this scene).  Contact forces then spike, FramesMO positions
        go to NaN, and the visual chain renders nothing -- the CTR
        appears to "disappear" right after the gap-profile PNG pops up.

        Fix: leave dt alone here.  The user changes dt from the Spinbox
        + Apply button at any time during the control phase, at their
        own pace.  default_control_dt is now used only as the Spinbox's
        initial *suggested* value, not as an applied value.
        """
        with self._lock:
            self._shared['init_complete'] = True
            self._shared['phase']         = 'control'

    def push_contact_profile(self, step, abscissae_m, gaps_m):
        """
        [ADDED] Called from SOFA thread by LiveContactMonitor every N steps.
        Lock-guarded write only -- NO Tk calls, NO matplotlib calls.

        The arrays are copied so that the caller's numpy buffers (which
        belong to BCM/contactMO and may be reallocated by SOFA at any time)
        cannot alias into our shared dict.

        Parameters
        ----------
        step : int
            Simulation step at which this frame was captured.  Used only
            for the plot title.
        abscissae_m : array-like of float
            Curvilinear abscissae of Pc_B along Tube_3 [meters], sorted.
            Empty array means "no valid contact pairs this frame" -- the
            plot will clear its markers.
        gaps_m : array-like of float
            Normal gap δn [meters] for each pair, in the same order as
            abscissae_m.  Negative values are interpenetration.
        """
        with self._lock:
            self._shared['contact_profile'] = {
                'step':      int(step),
                'abscissae': np.asarray(abscissae_m, dtype=float).copy(),
                'gaps_m':    np.asarray(gaps_m,      dtype=float).copy(),
                'dirty':     True,
            }

    # ==================================================================
    #  GUI thread entry point
    # ==================================================================
    def _run(self):
        try:
            self._tk_root = tk.Tk()
        except Exception as e:
            print(f"[CTRGui] FATAL: cannot create Tk root window: {e!r}")
            print( "[CTRGui] Simulation will run without GUI; use runSofa controls "
                   "and click Play in the runSofa toolbar.")
            return

        self._tk_root.title('CTR Two-Tube Controller')
        # Closing the window hides it without destroying; SOFA continues.
        self._tk_root.protocol("WM_DELETE_WINDOW",
                               lambda: self._tk_root.withdraw())

        self._build_widgets()

        # [ADDED] Live-plot Toplevel (separate window, child of main Tk root
        # so it shares the same event loop / mainloop / GIL state).  Two
        # tk.Tk() roots in one process do not share an event loop reliably;
        # a Toplevel does.
        if _MPL_OK:
            try:
                self._build_plot_panel()
                self._tk_root.after(self.PLOT_POLL_INTERVAL_MS, self._poll_plot)
            except Exception as e:
                print(f"[CTRGui] Live plot setup failed: {e!r}")

        self._tk_root.after(self.POLL_INTERVAL_MS, self._poll)

        try:
            self._tk_root.mainloop()
        except Exception as e:
            print(f"[CTRGui] Tk mainloop exited with exception: {e!r}")

    # ==================================================================
    #  Widget construction (all on GUI thread)
    # ==================================================================
    def _build_widgets(self):
        root = self._tk_root
        pad = 8

        # ---- Initialization frame ----
        top = ttk.LabelFrame(root, text='Initialization', padding=pad)
        top.pack(fill='x', padx=pad, pady=pad)

        self._init_button = ttk.Button(top, text='Initialize',
                                       command=self._on_init_clicked)
        self._init_button.pack(side='left', padx=(0, pad))

        self._status_label = ttk.Label(top, text='Status: waiting',
                                       font=('TkDefaultFont', 10, 'bold'))
        self._status_label.pack(side='left')

        # ---- External tube panel ----
        ext = ttk.LabelFrame(root, text='External tube  (Tube_1)', padding=pad)
        ext.pack(fill='x', padx=pad, pady=pad)
        self._control_frames.append(ext)
        (self._var_t1_tx, self._lbl_t1_tx,
         self._var_t1_rot, self._lbl_t1_rot) = self._build_tube_panel(
            ext,
            tx_max_cm=self._t1_max * 100.0,
            on_tx=self._on_t1_tx_change,
            on_rot=self._on_t1_rot_change,
            stop_cmd=lambda: self._reset_rot('t1'))

        # ---- Internal tube panel ----
        int_ = ttk.LabelFrame(root, text='Internal tube  (Tube_3)', padding=pad)
        int_.pack(fill='x', padx=pad, pady=pad)
        self._control_frames.append(int_)
        (self._var_t2_tx, self._lbl_t2_tx,
         self._var_t2_rot, self._lbl_t2_rot) = self._build_tube_panel(
            int_,
            tx_max_cm=self._t2_max * 100.0,
            on_tx=self._on_t2_tx_change,
            on_rot=self._on_t2_rot_change,
            stop_cmd=lambda: self._reset_rot('t2'))

        # ---- Time step & motion settings ----
        dt_frame = ttk.LabelFrame(root,
                                  text='Time step & motion  (control phase only)',
                                  padding=pad)
        dt_frame.pack(fill='x', padx=pad, pady=pad)
        self._control_frames.append(dt_frame)

        # ----- Row 0: dt -----
        ttk.Label(dt_frame, text='dt [s]:').grid(row=0, column=0, sticky='w')
        # The Spinbox starts at init_dt -- the actual dt the simulation is
        # running at the moment control phase begins.  Showing
        # default_control_dt here was misleading: the user could click Apply
        # without changing the value and trigger a silent 10x dt jump that
        # crashed the solver.
        self._var_dt = tk.StringVar(value=f"{self._init_dt:.6g}")
        sb = ttk.Spinbox(dt_frame, textvariable=self._var_dt, width=12,
                        from_=self._dt_min, to=self._dt_max,
                        increment=1e-4)
        sb.grid(row=0, column=1, padx=(pad, 0))
        ttk.Button(dt_frame, text='Apply',
                  command=self._on_dt_apply).grid(row=0, column=2, padx=(pad, 0))
        ttk.Label(dt_frame,
                 text=f'(suggested control dt = {self._default_control_dt:.0e};  '
                      f'allowed range [{self._dt_min:.0e}, {self._dt_max:.0e}])'
                 ).grid(row=0, column=3, padx=(pad, 0))

        # ----- Row 1: translation step size (live; no Apply button needed) -----
        # [MODIFIED, AGAIN] In *micrometers per simulation step*.  See the
        # docstring on shared['translation_step_m'] for why per-step rather
        # than per-second-sim-time is the right unit here.
        # Quick sanity check at default 50 µm/step:
        #   At ~100 sim steps per real second (typical for this scene),
        #   wall-clock speed is 50 µm * 100 = 5 mm / real second.
        #   To traverse 4 cm: ~8 real seconds.  To traverse 8 cm: ~16 s.
        # Bump up if you want faster, dial down if you want slower.
        ttk.Label(dt_frame, text='Translation step [µm/step]:'
                 ).grid(row=1, column=0, sticky='w', pady=(pad, 0))
        self._var_trans_step = tk.DoubleVar(value=self._default_trans_step_um)
        sp = ttk.Spinbox(dt_frame, textvariable=self._var_trans_step, width=12,
                        from_=self._trans_step_min_um,
                        to=self._trans_step_max_um,
                        increment=5.0,
                        command=self._on_trans_step_change)
        sp.grid(row=1, column=1, padx=(pad, 0), pady=(pad, 0))
        # Bind to keystrokes too, so typed-in values apply on Enter / focus-out.
        sp.bind('<Return>',   lambda _e: self._on_trans_step_change())
        sp.bind('<FocusOut>', lambda _e: self._on_trans_step_change())
        ttk.Label(dt_frame,
                 text=f'(default {self._default_trans_step_um:.0f} µm/step; '
                      f'range [{self._trans_step_min_um:.0f}, '
                      f'{self._trans_step_max_um:.0f}] µm/step; live)'
                 ).grid(row=1, column=2, columnspan=2,
                        sticky='w', padx=(pad, 0), pady=(pad, 0))

        # ----- Row 2: tip -----
        ttk.Label(dt_frame,
                 text='Tip: dt changes are now ramped (~2 %/step) toward the '
                      'target -- a 10x change takes ~120 steps, no more spasms.',
                 foreground='#666'
                 ).grid(row=2, column=0, columnspan=4,
                        sticky='w', pady=(4, 0))

        # Disable everything except the Initialize button until init completes.
        self._set_controls_state(False)

    def _build_tube_panel(self, parent, tx_max_cm, on_tx, on_rot, stop_cmd):
        """
        Tube panel layout (one frame per tube):

            Translation along X  [0 .. <max> cm]
            ╠════════════╪══════════════════════╣  <value> cm
            Rotation rate around X  [-<max> .. +<max> deg/step]   (+: CW, -: CCW)
            ╠══════════════╪══════════════════════╣  <value> deg/step
            [Stop rotation]
        """
        pad = 8

        ttk.Label(parent,
                 text=f'Translation along X  [0 .. {tx_max_cm:.0f} cm]'
                ).grid(row=0, column=0, sticky='w', pady=(0, 4))
        var_tx = tk.DoubleVar(value=0.0)
        ttk.Scale(parent, from_=0.0, to=tx_max_cm, orient='horizontal',
                 length=400, variable=var_tx, command=on_tx
                ).grid(row=1, column=0, sticky='ew')
        lbl_tx = ttk.Label(parent, text='0.00 cm', width=14, anchor='e')
        lbl_tx.grid(row=1, column=1, padx=(pad, 0))

        ttk.Label(parent,
                 text=f'Rotation rate around X  '
                      f'[{-self._max_rot_rate:+.1f} .. +{self._max_rot_rate:.1f} '
                      f'deg/step]   (+: CW,  -: CCW)'
                ).grid(row=2, column=0, sticky='w', pady=(pad, 4))
        var_rot = tk.DoubleVar(value=0.0)
        ttk.Scale(parent, from_=-self._max_rot_rate, to=self._max_rot_rate,
                 orient='horizontal', length=400, variable=var_rot,
                 command=on_rot
                ).grid(row=3, column=0, sticky='ew')
        lbl_rot = ttk.Label(parent, text='+0.00 deg/step', width=14, anchor='e')
        lbl_rot.grid(row=3, column=1, padx=(pad, 0))

        ttk.Button(parent, text='Stop rotation', command=stop_cmd
                  ).grid(row=4, column=0, sticky='w', pady=(pad, 0))

        return var_tx, lbl_tx, var_rot, lbl_rot

    # ==================================================================
    #  Enable / disable the control panel as a whole
    # ==================================================================
    def _set_controls_state(self, enabled):
        """
        Enable or disable every widget inside every control frame.
        ttk widgets use .state(); plain tk widgets use .configure(state=...).
        Some children (LabelFrame's internal Label, etc.) support neither
        and silently raise -- safe to ignore.
        """
        ttk_state    = ('!disabled',) if enabled else ('disabled',)
        plain_state  = 'normal'       if enabled else 'disabled'
        for fr in self._control_frames:
            for child in fr.winfo_children():
                try:
                    child.state(ttk_state)            # ttk widgets
                except (tk.TclError, AttributeError):
                    try:
                        child.configure(state=plain_state)  # plain tk widgets
                    except tk.TclError:
                        pass                                # unsupported; skip

    # ==================================================================
    #  Button / Scale handlers (run on GUI thread)
    # ==================================================================
    def _on_init_clicked(self):
        with self._lock:
            if self._shared['init_pressed']:
                return                                  # one-shot
            self._shared['init_pressed'] = True
            self._shared['phase']        = 'initializing'

        self._init_button.state(('disabled',))
        self._status_label.configure(text='Status: initializing (relaxing)…')

        # Unpause SOFA.  Single Data field write; benign if SOFA happens to
        # already have animate=True (e.g. user clicked runSofa Play earlier).
        try:
            self._root_node.animate = True
        except Exception as e:
            print(f"[CTRGui] Failed to unpause root_node.animate: {e!r}")

    def _on_t1_tx_change(self, _=None):
        v_cm = float(self._var_t1_tx.get())
        with self._lock:
            self._shared['t1_translation_target_m'] = v_cm * 0.01
        self._lbl_t1_tx.configure(text=f'{v_cm:.2f} cm')

    def _on_t1_rot_change(self, _=None):
        v = float(self._var_t1_rot.get())
        with self._lock:
            self._shared['t1_rotation_rate_deg_per_step'] = v
        self._lbl_t1_rot.configure(text=f'{v:+.2f} deg/step')

    def _on_t2_tx_change(self, _=None):
        v_cm = float(self._var_t2_tx.get())
        with self._lock:
            self._shared['t2_translation_target_m'] = v_cm * 0.01
        self._lbl_t2_tx.configure(text=f'{v_cm:.2f} cm')

    def _on_t2_rot_change(self, _=None):
        v = float(self._var_t2_rot.get())
        with self._lock:
            self._shared['t2_rotation_rate_deg_per_step'] = v
        self._lbl_t2_rot.configure(text=f'{v:+.2f} deg/step')

    def _reset_rot(self, which):
        if which == 't1':
            self._var_t1_rot.set(0.0)
            self._on_t1_rot_change()
        elif which == 't2':
            self._var_t2_rot.set(0.0)
            self._on_t2_rot_change()

    def _on_dt_apply(self):
        try:
            v = float(self._var_dt.get())
        except (ValueError, tk.TclError):
            print(f"[CTRGui] Invalid dt value: {self._var_dt.get()!r}")
            return
        if not (self._dt_min <= v <= self._dt_max):
            print(f"[CTRGui] dt={v:.6g} outside [{self._dt_min:.0e}, "
                  f"{self._dt_max:.0e}]; ignored.")
            return
        with self._lock:
            if self._shared['phase'] != 'control':
                print("[CTRGui] dt change requested but not in control phase; "
                      "ignored.")
                return
            self._shared['dt_request'] = v
        print(f"[CTRGui] dt change requested: {v:.6g} s "
              f"(applied at next AnimateBegin).")

    def _on_trans_step_change(self, _=None):
        """
        Live update of the translation per-step cap.  Called whenever the
        user clicks the Spinbox arrows, presses Enter, or tabs out of the
        field.  No Apply button -- the controller reads the shared value
        every step, so the new value takes effect on the next AnimateBegin.

        Units: GUI is µm/step; shared dict is m/step (controller's unit).

        Validation: clamp to [trans_step_min_um, trans_step_max_um] and
        reflect the clamped value back into the Spinbox so the user sees
        what was actually applied.  Bad / non-numeric input is rejected and
        the previous value is restored.
        """
        try:
            v_um = float(self._var_trans_step.get())
        except (ValueError, tk.TclError):
            # Restore last valid value
            with self._lock:
                v_m = float(self._shared['translation_step_m'])
            self._var_trans_step.set(v_m * 1e6)
            return

        clamped_um = max(self._trans_step_min_um,
                         min(self._trans_step_max_um, v_um))
        if clamped_um != v_um:
            self._var_trans_step.set(clamped_um)

        with self._lock:
            self._shared['translation_step_m'] = clamped_um * 1e-6

    # ==================================================================
    #  Polling loop  (runs on GUI thread, scheduled via root.after)
    # ==================================================================
    def _poll(self):
        """
        Runs every POLL_INTERVAL_MS on the GUI thread.

        Detects the init_complete flag set by the SOFA-thread call to
        signal_init_complete(); enables the control widgets exactly once,
        on the GUI thread, where Tk widget mutation is safe.
        """
        with self._lock:
            init_complete = self._shared['init_complete']

        if init_complete and not self._control_widgets_enabled:
            self._set_controls_state(True)
            self._status_label.configure(text='Status: control active')
            self._control_widgets_enabled = True

        # Reschedule -- Tk's after-callback model is the standard way to
        # implement periodic work without blocking the event loop.
        self._tk_root.after(self.POLL_INTERVAL_MS, self._poll)

    # ==================================================================
    #  Live contact-profile plot  (separate Toplevel window)
    # ==================================================================
    #
    # Why a Toplevel and not a second tk.Tk(): two Tk roots in one process
    # do not share an event loop reliably; the second mainloop tends to
    # starve the first.  A Toplevel is a child window of the main root --
    # one mainloop, one event-dispatch chain, two visible windows.
    #
    # Why FigureCanvasTkAgg instead of plt.show(): pyplot keeps a process-
    # global figure manager that fights with the embedding Tk; the
    # object-oriented Figure + FigureCanvasTkAgg pair has no global state
    # and behaves as a plain Tk widget.
    #
    # Why set_data / set_offsets and not re-plot: artist creation is the
    # slow path in matplotlib (allocating PathCollections, transforms,
    # etc.).  Updating the underlying numpy arrays of pre-built artists
    # and calling draw_idle() is what makes the plot feel "live" rather
    # than "stuttering at 2 fps".
    # ------------------------------------------------------------------
    def _build_plot_panel(self):
        """
        Create the live-plot Toplevel window with an embedded matplotlib
        canvas.  Called once on the GUI thread from _run().
        """
        win = tk.Toplevel(self._tk_root)
        win.title('CTR — Live contact profile')
        # Closing the window only hides it -- SOFA keeps stepping and the
        # buffer keeps being updated; the user can re-show with
        # self._plot_window.deiconify() if a re-show button is added later.
        win.protocol('WM_DELETE_WINDOW', win.withdraw)
        self._plot_window = win

        fig = Figure(figsize=(9.0, 4.5), dpi=100)
        ax = fig.add_subplot(111)

        # Reference line at δn = 0 (the wall-to-wall touching threshold).
        ax.axhline(0.0, color='k', linewidth=0.8, alpha=0.7)

        # Persistent artists -- their data is updated each frame by
        # _poll_plot(), no re-plot.
        (self._plot_line,) = ax.plot(
            [], [], '-', color='gray', linewidth=0.8, alpha=0.6)
        self._plot_scatter = ax.scatter(
            [], [], s=22, c=[], cmap='RdYlGn',
            vmin=-1.0, vmax=1.0,         # placeholder; rescaled per frame
            edgecolors='k', linewidths=0.4, zorder=3)

        ax.set_xlabel(r'Curvilinear abscissa of $P_{c,B}$ along Tube_3 [mm]')
        ax.set_ylabel(r'Normal gap $\delta_n$ [µm]')
        ax.set_title('Live contact profile  (waiting for data...)')
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        # Standard matplotlib pan/zoom/save toolbar -- nice for inspecting
        # an interesting frame without freezing the simulation.
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()
        canvas.draw()

        self._plot_canvas = canvas
        self._plot_ax     = ax

    def _poll_plot(self):
        """
        Plot poll loop.  Runs every PLOT_POLL_INTERVAL_MS on the Tk thread.
        Reads the latest frame from shared['contact_profile'] under lock,
        updates the matplotlib artists, requests a redraw.

        If the window has been hidden (user closed it via WM_DELETE), we
        keep polling -- the buffer is cheap and reopening should show
        fresh data.  If matplotlib threw at construction, this method was
        never scheduled in the first place, so we don't need to guard
        against missing artists here.
        """
        try:
            # Pull the latest frame.  Empty buffer / not-yet-written: skip.
            with self._lock:
                cp = self._shared['contact_profile']
                if not cp['dirty']:
                    return
                step      = cp['step']
                abscissae = cp['abscissae']
                gaps_m    = cp['gaps_m']
                cp['dirty'] = False

            # No-contact frame: clear all artists, keep axes labels.
            if abscissae is None or len(abscissae) == 0:
                self._plot_line.set_data([], [])
                self._plot_scatter.set_offsets(np.empty((0, 2)))
                self._plot_scatter.set_array(np.empty(0))
                self._plot_ax.set_title(
                    f'Live contact profile  (step {step}, K = 0)')
            else:
                # Convert SI -> human units for display.
                s_mm = abscissae * 1e3      # meters -> millimeters
                g_um = gaps_m    * 1e6      # meters -> micrometers

                # Update artist data in place.
                self._plot_line.set_data(s_mm, g_um)
                self._plot_scatter.set_offsets(np.c_[s_mm, g_um])
                self._plot_scatter.set_array(g_um)

                # Symmetric color-limit so green/red read as +/- gap sign.
                gmax = max(1e-12, float(np.abs(g_um).max()))
                self._plot_scatter.set_clim(-gmax, gmax)

                n_pen = int((g_um < 0).sum())
                self._plot_ax.set_title(
                    f'Live contact profile  (step {step}, '
                    f'K = {len(g_um)}, penetrating = {n_pen})')

                # Autoscale to current data each frame.  Cheap because
                # relim only walks the (small) artist list, not the data.
                self._plot_ax.relim()
                self._plot_ax.autoscale_view()

            # draw_idle is non-blocking and coalesces multiple requests --
            # the right primitive for live updates.
            self._plot_canvas.draw_idle()

        except Exception as e:
            # An exception inside an after-callback would silently kill
            # the periodic schedule.  Log + continue.
            print(f"[CTRGui] _poll_plot error: {e!r}")

        finally:
            # Reschedule unconditionally, even on error, so a transient
            # glitch does not permanently freeze the live plot.
            self._tk_root.after(self.PLOT_POLL_INTERVAL_MS, self._poll_plot)