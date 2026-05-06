# -*- coding: utf-8 -*-
"""
gui_5dof.py
===========
Tkinter-based GUI bridge for the curved-in-straight CTR scene
(ctr_curved_in_straight.py).

Differs from gui.py
-------------------
* The OUTER tube is fully pinned -- there are NO outer-tube controls.
* The INNER tube is driven from its proximal extremity with FIVE
  absolute-target sliders:

    tx [m]        translation along world +X        (range -max..+max)
    ty [m]        translation along world +Y        (range -max..+max)
    tz [m]        translation along world +Z        (range -max..+max)
    rx [deg]      rotation about local +X (TWIST)   (range -max..+max)
    rz [deg]      rotation about world  +Z (YAW)    (range -max..+max)

  Rotation about y is NOT exposed (the rigid base's PartialFixedProjective
  constraint pins it).

* Each axis is rate-limited per simulation step on the CONTROLLER side:
    translation_step_m   caps |tx,ty,tz| step size [m/step]
    rotation_step_rad    caps |rx,rz|    step size [rad/step]
  Both are read live from the snapshot every step.  Per-step (NOT
  per-sec-sim) is the right unit because the user perceives motion in
  real time, and real-time speed is roughly dt-independent for a
  non-adaptive solver.

Same external API as gui.CTRGuiBridge:
    snapshot()                -> dict copy of shared state (every step)
    consume_dt_request()      -> atomic read+clear of pending dt change
    signal_init_complete()    -> called by InitializationMonitor; flips
                                 phase to 'control', enables widgets
    push_contact_profile(...) -> called by LiveContactMonitor; updates
                                 the live-plot buffer

Threading architecture, phase state machine, live-plot buffering, and
the "GUI thread polls a flag, GUI thread mutates its own widgets"
pattern are IDENTICAL to gui.py.  See that file's module docstring for
the long-form rationale on Tk + SOFA + GIL coexistence.
"""

import math
import threading
import tkinter as tk
from tkinter import ttk

import numpy as np

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
    _MPL_OK = True
except ImportError as _e:
    print(f"[CTRGui5DOF] matplotlib not available, live plot disabled: {_e!r}")
    _MPL_OK = False


_NO_DT_REQUEST = None


class CTRGuiBridgeStraightOuter:
    """
    GUI bridge for the curved-in-straight CTR scene.

    Owns the Tkinter main window (5 sliders for the inner tube + dt /
    rate-step controls) and the live-plot Toplevel.  Construct in
    createScene() and stash on the CTRController so Python's GC can't
    collect it while SOFA runs.
    """

    POLL_INTERVAL_MS      = 100
    PLOT_POLL_INTERVAL_MS = 50

    # ==================================================================
    #  Construction
    # ==================================================================
    def __init__(self,
                 root_node,
                 # Inner-tube absolute-target ranges (sliders run -max..+max):
                 max_tx_m=0.04,
                 max_ty_m=0.005,
                 max_tz_m=0.005,
                 max_rx_deg=180.0,
                 max_rz_deg=15.0,
                 # Per-step rate limits + dt:
                 default_trans_step_um=50.0,
                 default_rot_step_deg=0.05,
                 init_dt=1e-4,
                 default_control_dt=1e-3,
                 dt_min=1e-6, dt_max=1e-1,
                 trans_step_min_um=1.0,
                 trans_step_max_um=500.0,
                 rot_step_min_deg=0.001,
                 rot_step_max_deg=2.0):

        self._root_node            = root_node
        self._max_tx_m             = float(max_tx_m)
        self._max_ty_m             = float(max_ty_m)
        self._max_tz_m             = float(max_tz_m)
        self._max_rx_deg           = float(max_rx_deg)
        self._max_rz_deg           = float(max_rz_deg)
        self._init_dt              = float(init_dt)
        self._default_control_dt   = float(default_control_dt)
        self._dt_min               = float(dt_min)
        self._dt_max               = float(dt_max)
        self._default_trans_step_um = float(default_trans_step_um)
        self._default_rot_step_deg  = float(default_rot_step_deg)
        self._trans_step_min_um    = float(trans_step_min_um)
        self._trans_step_max_um    = float(trans_step_max_um)
        self._rot_step_min_deg     = float(rot_step_min_deg)
        self._rot_step_max_deg     = float(rot_step_max_deg)

        # ---- Shared state (guarded by self._lock; readable from any thread) ----
        self._lock = threading.Lock()
        self._shared = {
            'phase':                  'waiting',
            'init_pressed':           False,
            'init_complete':          False,
            # Inner-tube absolute targets (5 DOFs).  All in SI radians/meters
            # internally; the GUI displays cm/deg but writes m/rad here.
            't2_tx_target_m':         0.0,
            't2_ty_target_m':         0.0,
            't2_tz_target_m':         0.0,
            't2_rx_target_rad':       0.0,
            't2_rz_target_rad':       0.0,
            # Per-step rate caps -- read live by the controller every step.
            'translation_step_m':     float(default_trans_step_um) * 1e-6,
            'rotation_step_rad':      math.radians(float(default_rot_step_deg)),
            'dt_request':             _NO_DT_REQUEST,
            # Live contact-profile buffer (single-frame "latest value" slot).
            # Same semantics as gui.py: SOFA writes via push_contact_profile,
            # GUI reads via _poll_plot, dirty flag clears on consume.
            'contact_profile': {
                'step':      0,
                'abscissae': None,
                'gaps_m':    None,
                'dirty':     False,
            },
        }

        # ---- GUI-thread-only state ----
        self._tk_root                 = None
        self._init_button             = None
        self._status_label            = None
        self._control_frames          = []
        self._control_widgets_enabled = False

        # Tk variables (one per slider/spinbox)
        self._var_tx = self._var_ty = self._var_tz = None
        self._var_rx = self._var_rz = None
        self._var_dt          = None
        self._var_trans_step  = None        # micrometers/step
        self._var_rot_step    = None        # degrees/step

        # Live value labels
        self._lbl_tx = self._lbl_ty = self._lbl_tz = None
        self._lbl_rx = self._lbl_rz = None

        # Live-plot Tk Toplevel + matplotlib state (GUI-thread only)
        self._plot_window  = None
        self._plot_canvas  = None
        self._plot_ax      = None
        self._plot_line    = None
        self._plot_scatter = None

        # Start the GUI thread.  Daemon -> dies with the SOFA process.
        self._thread = threading.Thread(target=self._run,
                                        name='CTRGui5DOFThread',
                                        daemon=True)
        self._thread.start()

    # ==================================================================
    #  SOFA-side API  (thread-safe; callable from any thread)
    # ==================================================================
    def snapshot(self):
        """Atomic dict copy of the shared state.  Called every step."""
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
        Flips phase to 'control' and sets the init_complete flag the
        GUI poll loop is waiting on.  dt is NOT changed here -- the
        ctr_two_tubes.py post-mortem (auto-bumping dt at this point
        crashed the constraint solve) applies equally here.  See
        gui.py:signal_init_complete for the full rationale.
        """
        with self._lock:
            self._shared['phase']         = 'control'
            self._shared['init_complete'] = True

    def push_contact_profile(self, step, abscissae_m, gaps_m):
        """
        Called by LiveContactMonitor every N simulation steps.  Replaces
        the single-frame 'latest value' slot in the shared dict.  No
        queue -- if the Tk repaint is slower than the SOFA push rate,
        intermediate frames are silently overwritten (correct behaviour
        for a live monitor: always show the freshest data).
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
            print(f"[CTRGui5DOF] FATAL: cannot create Tk root window: {e!r}")
            print( "[CTRGui5DOF] Simulation will run without GUI; use runSofa "
                   "controls and click Play in the runSofa toolbar.")
            return

        self._tk_root.title('CTR Curved-in-Straight Controller')
        # Closing the window only hides it; SOFA continues.
        self._tk_root.protocol("WM_DELETE_WINDOW",
                               lambda: self._tk_root.withdraw())

        self._build_widgets()

        if _MPL_OK:
            try:
                self._build_plot_panel()
                self._tk_root.after(self.PLOT_POLL_INTERVAL_MS, self._poll_plot)
            except Exception as e:
                print(f"[CTRGui5DOF] Live plot setup failed: {e!r}")

        self._tk_root.after(self.POLL_INTERVAL_MS, self._poll)

        try:
            self._tk_root.mainloop()
        except Exception as e:
            print(f"[CTRGui5DOF] Tk mainloop exited with exception: {e!r}")

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

        # ---- Inner-tube panel (5 absolute-target sliders) ----
        inner = ttk.LabelFrame(root,
                               text='Inner tube  (proximal extremity, '
                                    '5-DOF absolute targets)',
                               padding=pad)
        inner.pack(fill='x', padx=pad, pady=pad)
        self._control_frames.append(inner)

        # Translation block
        tx_max_cm = self._max_tx_m * 100.0
        ty_max_cm = self._max_ty_m * 100.0
        tz_max_cm = self._max_tz_m * 100.0

        (self._var_tx, self._lbl_tx) = self._build_translation_slider(
            inner, row=0, axis_name='X',
            max_cm=tx_max_cm,
            on_change=self._on_tx_change)
        (self._var_ty, self._lbl_ty) = self._build_translation_slider(
            inner, row=2, axis_name='Y',
            max_cm=ty_max_cm,
            on_change=self._on_ty_change)
        (self._var_tz, self._lbl_tz) = self._build_translation_slider(
            inner, row=4, axis_name='Z',
            max_cm=tz_max_cm,
            on_change=self._on_tz_change)

        # Rotation block
        (self._var_rx, self._lbl_rx) = self._build_rotation_slider(
            inner, row=6, axis_name='X (twist)',
            max_deg=self._max_rx_deg,
            on_change=self._on_rx_change)
        (self._var_rz, self._lbl_rz) = self._build_rotation_slider(
            inner, row=8, axis_name='Z (yaw)',
            max_deg=self._max_rz_deg,
            on_change=self._on_rz_change)

        # Reset-all button (handy for snapping back to home pose)
        ttk.Button(inner, text='Reset all DOFs to 0',
                   command=self._on_reset_all
                  ).grid(row=10, column=0, sticky='w', pady=(pad, 0))

        # ---- Time step & motion settings ----
        dt_frame = ttk.LabelFrame(root,
                                  text='Time step & motion  (control phase only)',
                                  padding=pad)
        dt_frame.pack(fill='x', padx=pad, pady=pad)
        self._control_frames.append(dt_frame)

        # Row 0: dt
        ttk.Label(dt_frame, text='dt [s]:').grid(row=0, column=0, sticky='w')
        self._var_dt = tk.StringVar(value=f"{self._init_dt:.6g}")
        ttk.Spinbox(dt_frame, textvariable=self._var_dt, width=12,
                    from_=self._dt_min, to=self._dt_max,
                    increment=1e-4
                   ).grid(row=0, column=1, padx=(pad, 0))
        ttk.Button(dt_frame, text='Apply',
                   command=self._on_dt_apply
                  ).grid(row=0, column=2, padx=(pad, 0))
        ttk.Label(dt_frame,
                  text=f'(suggested control dt = '
                       f'{self._default_control_dt:.0e};  '
                       f'range [{self._dt_min:.0e}, {self._dt_max:.0e}])'
                 ).grid(row=0, column=3, padx=(pad, 0))

        # Row 1: translation step
        ttk.Label(dt_frame, text='Translation step [µm/step]:'
                 ).grid(row=1, column=0, sticky='w', pady=(pad, 0))
        self._var_trans_step = tk.DoubleVar(value=self._default_trans_step_um)
        sp = ttk.Spinbox(dt_frame, textvariable=self._var_trans_step, width=12,
                         from_=self._trans_step_min_um,
                         to=self._trans_step_max_um,
                         increment=5.0,
                         command=self._on_trans_step_change)
        sp.grid(row=1, column=1, padx=(pad, 0), pady=(pad, 0))
        sp.bind('<Return>',   lambda _e: self._on_trans_step_change())
        sp.bind('<FocusOut>', lambda _e: self._on_trans_step_change())
        ttk.Label(dt_frame,
                  text=f'(default {self._default_trans_step_um:.0f} µm/step; '
                       f'range [{self._trans_step_min_um:.0f}, '
                       f'{self._trans_step_max_um:.0f}]; live)'
                 ).grid(row=1, column=2, columnspan=2,
                        sticky='w', padx=(pad, 0), pady=(pad, 0))

        # Row 2: rotation step
        ttk.Label(dt_frame, text='Rotation step [deg/step]:'
                 ).grid(row=2, column=0, sticky='w', pady=(pad, 0))
        self._var_rot_step = tk.DoubleVar(value=self._default_rot_step_deg)
        sp_r = ttk.Spinbox(dt_frame, textvariable=self._var_rot_step, width=12,
                           from_=self._rot_step_min_deg,
                           to=self._rot_step_max_deg,
                           increment=0.01,
                           command=self._on_rot_step_change)
        sp_r.grid(row=2, column=1, padx=(pad, 0), pady=(pad, 0))
        sp_r.bind('<Return>',   lambda _e: self._on_rot_step_change())
        sp_r.bind('<FocusOut>', lambda _e: self._on_rot_step_change())
        ttk.Label(dt_frame,
                  text=f'(default {self._default_rot_step_deg:.3f} deg/step; '
                       f'range [{self._rot_step_min_deg:.3f}, '
                       f'{self._rot_step_max_deg:.3f}]; live)'
                 ).grid(row=2, column=2, columnspan=2,
                        sticky='w', padx=(pad, 0), pady=(pad, 0))

        # Row 3: dt-ramp tip
        ttk.Label(dt_frame,
                  text='Tip: dt changes are ramped (~2 %/step) toward the '
                       'target; per-step caps make motion dt-independent.',
                  foreground='#666'
                 ).grid(row=3, column=0, columnspan=4,
                        sticky='w', pady=(4, 0))

        # Disable everything except the Initialize button until init completes.
        self._set_controls_state(False)

    def _build_translation_slider(self, parent, row, axis_name, max_cm, on_change):
        """One translation slider with -max_cm..+max_cm range."""
        pad = 8
        ttk.Label(parent,
                  text=f'Translation along {axis_name}  '
                       f'[{-max_cm:+.2f} .. +{max_cm:.2f} cm]'
                 ).grid(row=row, column=0, sticky='w', pady=(0, 4))
        var = tk.DoubleVar(value=0.0)
        ttk.Scale(parent, from_=-max_cm, to=max_cm, orient='horizontal',
                  length=400, variable=var, command=on_change
                 ).grid(row=row + 1, column=0, sticky='ew')
        lbl = ttk.Label(parent, text='+0.00 cm', width=14, anchor='e')
        lbl.grid(row=row + 1, column=1, padx=(pad, 0))
        return var, lbl

    def _build_rotation_slider(self, parent, row, axis_name, max_deg, on_change):
        """One rotation slider with -max_deg..+max_deg range."""
        pad = 8
        ttk.Label(parent,
                  text=f'Rotation about {axis_name}  '
                       f'[{-max_deg:+.1f} .. +{max_deg:.1f} deg]'
                 ).grid(row=row, column=0, sticky='w', pady=(pad, 4))
        var = tk.DoubleVar(value=0.0)
        ttk.Scale(parent, from_=-max_deg, to=max_deg, orient='horizontal',
                  length=400, variable=var, command=on_change
                 ).grid(row=row + 1, column=0, sticky='ew')
        lbl = ttk.Label(parent, text='+0.00 deg', width=14, anchor='e')
        lbl.grid(row=row + 1, column=1, padx=(pad, 0))
        return var, lbl

    # ==================================================================
    #  Enable / disable the control panel as a whole
    # ==================================================================
    def _set_controls_state(self, enabled):
        """Identical to gui.py: ttk widgets use .state(); plain tk widgets
        use .configure(state=...).  Some children silently raise; ignore."""
        ttk_state    = ('!disabled',) if enabled else ('disabled',)
        plain_state  = 'normal'       if enabled else 'disabled'
        for fr in self._control_frames:
            for child in fr.winfo_children():
                try:
                    child.state(ttk_state)
                except (tk.TclError, AttributeError):
                    try:
                        child.configure(state=plain_state)
                    except tk.TclError:
                        pass

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

        try:
            self._root_node.animate = True
        except Exception as e:
            print(f"[CTRGui5DOF] Failed to unpause root_node.animate: {e!r}")

    def _on_tx_change(self, _=None):
        v_cm = float(self._var_tx.get())
        with self._lock:
            self._shared['t2_tx_target_m'] = v_cm * 0.01
        self._lbl_tx.configure(text=f'{v_cm:+.2f} cm')

    def _on_ty_change(self, _=None):
        v_cm = float(self._var_ty.get())
        with self._lock:
            self._shared['t2_ty_target_m'] = v_cm * 0.01
        self._lbl_ty.configure(text=f'{v_cm:+.2f} cm')

    def _on_tz_change(self, _=None):
        v_cm = float(self._var_tz.get())
        with self._lock:
            self._shared['t2_tz_target_m'] = v_cm * 0.01
        self._lbl_tz.configure(text=f'{v_cm:+.2f} cm')

    def _on_rx_change(self, _=None):
        v_deg = float(self._var_rx.get())
        with self._lock:
            self._shared['t2_rx_target_rad'] = math.radians(v_deg)
        self._lbl_rx.configure(text=f'{v_deg:+.2f} deg')

    def _on_rz_change(self, _=None):
        v_deg = float(self._var_rz.get())
        with self._lock:
            self._shared['t2_rz_target_rad'] = math.radians(v_deg)
        self._lbl_rz.configure(text=f'{v_deg:+.2f} deg')

    def _on_reset_all(self):
        """Snap all 5 sliders + their shared targets back to 0."""
        for var, on_change in (
            (self._var_tx, self._on_tx_change),
            (self._var_ty, self._on_ty_change),
            (self._var_tz, self._on_tz_change),
            (self._var_rx, self._on_rx_change),
            (self._var_rz, self._on_rz_change),
        ):
            var.set(0.0)
            on_change()

    def _on_dt_apply(self):
        try:
            v = float(self._var_dt.get())
        except (ValueError, tk.TclError):
            print(f"[CTRGui5DOF] Invalid dt value: {self._var_dt.get()!r}")
            return
        if not (self._dt_min <= v <= self._dt_max):
            print(f"[CTRGui5DOF] dt={v:.6g} outside "
                  f"[{self._dt_min:.0e}, {self._dt_max:.0e}]; ignored.")
            return
        with self._lock:
            if self._shared['phase'] != 'control':
                print("[CTRGui5DOF] dt change requested but not in control "
                      "phase; ignored.")
                return
            self._shared['dt_request'] = v
        print(f"[CTRGui5DOF] dt request queued: {v:.6g} s "
              f"(controller will ramp toward this).")

    def _on_trans_step_change(self):
        try:
            v_um = float(self._var_trans_step.get())
        except (ValueError, tk.TclError):
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

    def _on_rot_step_change(self):
        try:
            v_deg = float(self._var_rot_step.get())
        except (ValueError, tk.TclError):
            with self._lock:
                v_rad = float(self._shared['rotation_step_rad'])
            self._var_rot_step.set(math.degrees(v_rad))
            return

        clamped_deg = max(self._rot_step_min_deg,
                          min(self._rot_step_max_deg, v_deg))
        if clamped_deg != v_deg:
            self._var_rot_step.set(clamped_deg)

        with self._lock:
            self._shared['rotation_step_rad'] = math.radians(clamped_deg)

    # ==================================================================
    #  Polling loop  (runs on GUI thread, scheduled via root.after)
    # ==================================================================
    def _poll(self):
        """Same pattern as gui.py: detect init_complete flag set by the
        SOFA-thread call to signal_init_complete, enable the control
        widgets exactly once on the GUI thread (where Tk widget mutation
        is safe)."""
        with self._lock:
            init_complete = self._shared['init_complete']

        if init_complete and not self._control_widgets_enabled:
            self._set_controls_state(True)
            self._status_label.configure(text='Status: control active')
            self._control_widgets_enabled = True

        self._tk_root.after(self.POLL_INTERVAL_MS, self._poll)

    # ==================================================================
    #  Live contact-profile plot  (separate Toplevel window)
    # ==================================================================
    # See gui.py:_build_plot_panel for the long-form design rationale.
    # ------------------------------------------------------------------
    def _build_plot_panel(self):
        win = tk.Toplevel(self._tk_root)
        win.title('CTR — Live contact profile')
        win.protocol('WM_DELETE_WINDOW', win.withdraw)
        self._plot_window = win

        fig = Figure(figsize=(9.0, 4.5), dpi=100)
        ax = fig.add_subplot(111)
        ax.axhline(0.0, color='k', linewidth=0.8, alpha=0.7)

        (self._plot_line,) = ax.plot(
            [], [], '-', color='gray', linewidth=0.8, alpha=0.6)
        self._plot_scatter = ax.scatter(
            [], [], s=22, c=[], cmap='RdYlGn',
            vmin=-1.0, vmax=1.0,
            edgecolors='k', linewidths=0.4, zorder=3)

        ax.set_xlabel(r'Curvilinear abscissa of $P_{c,B}$ along inner tube [mm]')
        ax.set_ylabel(r'Normal gap $\delta_n$ [µm]')
        ax.set_title('Live contact profile  (waiting for data...)')
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()
        canvas.draw()

        self._plot_canvas = canvas
        self._plot_ax     = ax

    def _poll_plot(self):
        """Same logic as gui.py:_poll_plot -- pull latest frame under
        lock, update artists in place, redraw via draw_idle."""
        try:
            with self._lock:
                cp = self._shared['contact_profile']
                if not cp['dirty']:
                    return
                step      = cp['step']
                abscissae = cp['abscissae']
                gaps_m    = cp['gaps_m']
                cp['dirty'] = False

            if abscissae is None or len(abscissae) == 0:
                self._plot_line.set_data([], [])
                self._plot_scatter.set_offsets(np.empty((0, 2)))
                self._plot_scatter.set_array(np.empty(0))
                self._plot_ax.set_title(
                    f'Live contact profile  (step {step}, K = 0)')
            else:
                s_mm = abscissae * 1e3
                g_um = gaps_m    * 1e6

                self._plot_line.set_data(s_mm, g_um)
                self._plot_scatter.set_offsets(np.c_[s_mm, g_um])
                self._plot_scatter.set_array(g_um)

                gmax = max(1e-12, float(np.abs(g_um).max()))
                self._plot_scatter.set_clim(-gmax, gmax)

                n_pen = int((g_um < 0).sum())
                self._plot_ax.set_title(
                    f'Live contact profile  (step {step}, '
                    f'K = {len(g_um)}, penetrating = {n_pen})')

                self._plot_ax.relim()
                self._plot_ax.autoscale_view()

            self._plot_canvas.draw_idle()

        except Exception as e:
            print(f"[CTRGui5DOF] _poll_plot error: {e!r}")

        finally:
            self._tk_root.after(self.PLOT_POLL_INTERVAL_MS, self._poll_plot)