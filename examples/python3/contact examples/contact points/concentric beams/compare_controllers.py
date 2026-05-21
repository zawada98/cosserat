# -*- coding: utf-8 -*-
"""
compare_controllers.py
================================================================================
Headless-comparison controllers and schedule helper for benchmarking two
constraint-class variants of the two-tube CTR scene (CPULC vs stock ULC + feeder).

Three classes:

    Schedule              -- pure-Python state machine over the simulation
                              step counter; encodes init -> translation
                              increments -> rotation increments -> done.
                              Used by both the actuator (decide what pose to
                              write) and the recorder (decide when to take a
                              heavy snapshot).

    ScriptedActuator      -- Sofa.Core.Controller; replaces the GUI-driven
                              CTRController.  Reads the schedule each step and
                              writes the rigid base poses of Tube_1 (always
                              zero) and Tube_3 (cumulative) using the static
                              method CTRController._set_pose imported from
                              ctr_two_tubes.

    ComparisonRecorder    -- Sofa.Core.Controller; replaces LiveContactMonitor.
                              Captures per-step diagnostics (every step) and
                              per-snapshot heavy data (gap/abscissa profile)
                              at end-of-hold steps and end-of-init.  Writes a
                              single .npz at end of run (also via atexit, in
                              case the run is killed early).

The two controller classes are pure SOFA-thread; no Tk, no matplotlib, no
threads.  All output is to a .npz file.

The contact-MO read pattern (Pc_A/Pc_B interleaved at index 2k+1) and the
polyline projection of Pc_B onto Tube_3's centerline are lifted verbatim from
LiveContactMonitor (live_monitor.py).
"""

from __future__ import annotations

import atexit
import math
import os
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
import Sofa
import Sofa.Core


# Sentinel from BeamContactMapping.cpp (kInvalidGap = 1e9) marks unused slots.
INVALID_GAP_THRESHOLD = 1e8


# =============================================================================
#  Schedule
# =============================================================================

@dataclass
class Schedule:
    """
    Pure-Python state machine over the simulation step counter.

    Phases (in order):
        'init'       0 .. INIT_STEPS-1               inner tube relaxes
        'trans'      next TRANS_N_INCREMENTS *       30 (1 mm + hold) cycles
                       (TRANS_STEPS_PER_INCR + HOLD_STEPS) steps
        'rot'        next ROT_N_INCREMENTS *         30 (1 deg + hold) cycles
                       (ROT_STEPS_PER_INCR + HOLD_STEPS) steps
        'done'       any step >= total_steps         past the end of run

    Each translation/rotation cycle is split into 'actuating' and 'holding'
    sub-phases.  A snapshot is fired at the LAST step of each holding phase
    (and at the last step of init), giving 1 + 30 + 30 = 61 snapshot points.
    """
    init_steps:           int
    hold_steps:           int
    trans_n_increments:   int
    trans_steps_per_incr: int
    rot_n_increments:     int
    rot_steps_per_incr:   int

    # ------------------------------------------------------------------
    @property
    def trans_cycle_len(self) -> int:
        return self.trans_steps_per_incr + self.hold_steps

    @property
    def rot_cycle_len(self) -> int:
        return self.rot_steps_per_incr + self.hold_steps

    @property
    def trans_phase_start(self) -> int:
        return self.init_steps

    @property
    def rot_phase_start(self) -> int:
        return self.trans_phase_start + self.trans_n_increments * self.trans_cycle_len

    @property
    def total_steps(self) -> int:
        return self.rot_phase_start + self.rot_n_increments * self.rot_cycle_len

    @property
    def total_snapshots(self) -> int:
        return 1 + self.trans_n_increments + self.rot_n_increments

    # ------------------------------------------------------------------
    def phase_of(self, step: int) -> dict:
        """
        Resolve `step` (1-based; first sim step is step=1) into a phase dict:
            {phase: 'init'|'trans'|'rot'|'done',
             sub:   'actuating'|'holding'|'-',
             incr:  int,           # 0..N-1 within trans/rot, else -1
             step_in_incr: int}    # 0-based position within the current cycle

        Step counter convention: ScriptedActuator increments self._step at
        the start of onAnimateBeginEvent BEFORE calling this, so the first
        physics step is step=1.  Phase boundaries below are written so that
        the last step of the init phase is step=init_steps.
        """
        if step <= 0:
            return {'phase': 'init', 'sub': '-', 'incr': -1, 'step_in_incr': 0}

        if step <= self.init_steps:
            return {'phase': 'init', 'sub': '-', 'incr': -1,
                    'step_in_incr': step - 1}

        # In trans phase?
        s = step - self.trans_phase_start                 # 1..N*cycle_len
        if s <= self.trans_n_increments * self.trans_cycle_len:
            incr = (s - 1) // self.trans_cycle_len        # 0..N-1
            within = (s - 1) % self.trans_cycle_len       # 0..cycle_len-1
            sub = 'actuating' if within < self.trans_steps_per_incr else 'holding'
            return {'phase': 'trans', 'sub': sub, 'incr': incr,
                    'step_in_incr': within}

        # In rot phase?
        s = step - self.rot_phase_start
        if s <= self.rot_n_increments * self.rot_cycle_len:
            incr = (s - 1) // self.rot_cycle_len
            within = (s - 1) % self.rot_cycle_len
            sub = 'actuating' if within < self.rot_steps_per_incr else 'holding'
            return {'phase': 'rot', 'sub': sub, 'incr': incr,
                    'step_in_incr': within}

        return {'phase': 'done', 'sub': '-', 'incr': -1, 'step_in_incr': 0}

    # ------------------------------------------------------------------
    def is_snapshot_step(self, step: int) -> bool:
        """
        True iff `step` is the last step of either the init phase or any
        hold phase.  This is the moment to capture a heavy snapshot.
        """
        if step == self.init_steps:
            return True
        # Last step of any trans-cycle hold:
        if self.trans_phase_start < step <= self.rot_phase_start:
            offset = step - self.trans_phase_start
            return (offset % self.trans_cycle_len) == 0
        # Last step of any rot-cycle hold:
        if self.rot_phase_start < step <= self.total_steps:
            offset = step - self.rot_phase_start
            return (offset % self.rot_cycle_len) == 0
        return False

    # ------------------------------------------------------------------
    def snapshot_index(self, step: int) -> int:
        """
        Map a snapshot step to its index in [0 .. total_snapshots-1].
        Returns -1 if `step` is not a snapshot step.

        Layout:
            0                                  = end-of-init
            1 .. N_trans                       = end of each trans hold
            N_trans+1 .. N_trans+N_rot         = end of each rot hold
        """
        if not self.is_snapshot_step(step):
            return -1
        if step == self.init_steps:
            return 0
        if step <= self.rot_phase_start:
            return (step - self.trans_phase_start) // self.trans_cycle_len
        return (self.trans_n_increments
                + (step - self.rot_phase_start) // self.rot_cycle_len)


# =============================================================================
#  ScriptedActuator
# =============================================================================

class ScriptedActuator(Sofa.Core.Controller):
    """
    Drop-in replacement for CTRController in the headless comparison run.

    No GUI, no phase machine reading from a bridge.  Each step it consults
    the Schedule and integrates the cumulative Tube_3 base pose accordingly:

      'init', 'holding'   -> no change to accumulators
      'trans/actuating'   -> _t2_pos_m   += trans_delta_per_step
      'rot/actuating'     -> _t2_angle  += rot_delta_per_step

    Tube_1 base is held at (0, 0) throughout.  Both bases are written every
    step via the same _set_pose helper used by CTRController, which we
    import from ctr_two_tubes.
    """

    def __init__(self,
                 root_node,
                 t1_base_mo,
                 t2_base_mo,
                 t2_x_offset,
                 schedule: Schedule,
                 trans_delta_per_step: float,    # m / step
                 rot_delta_per_step:   float,    # rad / step
                 set_pose_fn,                    # CTRController._set_pose
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.root_node = root_node
        self.t1_base_mo = t1_base_mo
        self.t2_base_mo = t2_base_mo
        self.t2_x0 = float(t2_x_offset)
        self.schedule = schedule
        self.dx = float(trans_delta_per_step)
        self.dth = float(rot_delta_per_step)
        self._set_pose = set_pose_fn

        # Cumulative Tube_3 pose, never reset.
        self._t2_pos_m = 0.0
        self._t2_angle_rad = 0.0

        # Step counter: pre-incremented at start of begin-event, so first
        # physics step is step=1 (matches schedule.phase_of convention).
        self._step = 0

        # Cached previous phase string for end-of-phase logging.
        self._prev_phase = 'init'

    # ------------------------------------------------------------------
    def onAnimateBeginEvent(self, event):
        self._step += 1
        ph = self.schedule.phase_of(self._step)

        # Phase-transition log lines, useful when tailing the runSofa stdout.
        phase_label = f"{ph['phase']}/{ph['sub']}"
        if phase_label != self._prev_phase:
            print(f"[ScriptedActuator] step={self._step} "
                  f"phase {self._prev_phase} -> {phase_label}  "
                  f"t2_pos={self._t2_pos_m*1e3:.3f} mm  "
                  f"t2_angle={math.degrees(self._t2_angle_rad):.3f} deg",
                  flush=True)
            self._prev_phase = phase_label

        # Integrate accumulators only during actuation sub-phases.
        if ph['phase'] == 'trans' and ph['sub'] == 'actuating':
            self._t2_pos_m += self.dx
        elif ph['phase'] == 'rot' and ph['sub'] == 'actuating':
            self._t2_angle_rad += self.dth
        # 'init', 'holding', 'done' -> hold accumulators

        # Tube_1 base: zero displacement, zero rotation, every step.
        self._set_pose(self.t1_base_mo, 0.0, 0.0, x0=0.0)
        # Tube_3 base: cumulative pose with intrinsic offset.
        self._set_pose(self.t2_base_mo, self._t2_pos_m, self._t2_angle_rad,
                       x0=self.t2_x0)

        # Past end of schedule: stop animating (the recorder writes the
        # final npz on its own end-event for this step; runSofa --start in
        # batch mode then idles until the driver process kills it).
        if ph['phase'] == 'done':
            try:
                self.root_node.animate = False
            except Exception:
                pass


# =============================================================================
#  ComparisonRecorder
# =============================================================================

class ComparisonRecorder(Sofa.Core.Controller):
    """
    Captures comparison data and writes a single .npz at end of run.

    Per-step (every step):
        step_idx, sim_time, t2_pos_m, t2_angle_rad,
        n_iters, residual, max_pen, n_active, wall_dt_s

    Per-snapshot (61 entries, at end-of-init and end-of-each-hold):
        snap_step_idx, snap_phase_label, snap_t2_pos_m, snap_t2_angle_rad,
        snap_s_B (variable-length float array per snap),
        snap_gap_n (variable-length float array per snap)

    Solver field probing
    --------------------
    BlockGaussSeidelConstraintSolver may expose its current GS iteration
    count and residual under any of several names depending on plugin
    version.  At construction time we probe a list of candidate Data field
    names; the first that resolves is used.  If none resolve, NaN is
    recorded and a one-time warning lists what fields ARE available so the
    user can correct the candidate list.
    """

    SOLVER_ITER_CANDIDATES = [
        'currentIterations', 'currentNumIterations',
        'iterations', 'iter', 'numIter', 'd_currentIterations',
    ]
    SOLVER_ERROR_CANDIDATES = [
        'currentError', 'currentErrorNorm',
        'error', 'residual', 'd_currentError',
    ]

    def __init__(self,
                 root_node,
                 actuator: ScriptedActuator,
                 schedule: Schedule,
                 constraint_solver,
                 t2_MO,
                 bcm,
                 contact_mo,
                 t2_frame_curv_abs,
                 out_path: str,
                 mode_label: str,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.root_node = root_node
        self.actuator = actuator
        self.schedule = schedule
        self.solver = constraint_solver
        self.t2_MO = t2_MO
        self.bcm = bcm
        self.contact_mo = contact_mo
        self.t2_curv_abs = np.asarray(t2_frame_curv_abs, dtype=float)
        self.out_path = out_path
        self.mode_label = mode_label

        # ---- Probe solver fields once at construction --------------------
        self._iter_field = self._probe_field(self.solver, self.SOLVER_ITER_CANDIDATES)
        self._error_field = self._probe_field(self.solver, self.SOLVER_ERROR_CANDIDATES)
        if self._iter_field is None or self._error_field is None:
            avail = self._list_available_data(self.solver)
            print(f"[ComparisonRecorder] WARNING: solver probe missed "
                  f"(iter={self._iter_field}, error={self._error_field}). "
                  f"Available Data fields on {type(self.solver).__name__}:\n"
                  f"  {avail}", flush=True)
        else:
            print(f"[ComparisonRecorder] solver fields resolved: "
                  f"iter='{self._iter_field}', error='{self._error_field}'",
                  flush=True)

        # ---- Pre-allocate per-step buffers up to total_steps + small margin
        T = schedule.total_steps + 16
        self._step_idx     = np.full(T, -1,  dtype=np.int64)
        self._sim_time     = np.full(T, np.nan, dtype=np.float64)
        self._t2_pos_m     = np.full(T, np.nan, dtype=np.float64)
        self._t2_angle_rad = np.full(T, np.nan, dtype=np.float64)
        self._n_iters      = np.full(T, np.nan, dtype=np.float64)
        self._residual     = np.full(T, np.nan, dtype=np.float64)
        self._max_pen      = np.full(T, np.nan, dtype=np.float64)
        self._n_active     = np.full(T, -1,  dtype=np.int64)
        self._wall_dt_s    = np.full(T, np.nan, dtype=np.float64)

        # Per-snapshot lists (variable-length s/gap arrays).
        self._snap_step    = []
        self._snap_phase   = []
        self._snap_t2_pos  = []
        self._snap_t2_ang  = []
        self._snap_s_B     = []
        self._snap_gap_n   = []

        self._step  = 0
        self._t_begin = None        # set in onAnimateBeginEvent
        self._written = False

        # Backup write on Python exit (in case runSofa is killed).
        atexit.register(self._write_npz)

    # ------------------------------------------------------------------
    def onAnimateBeginEvent(self, event):
        self._t_begin = time.perf_counter()

    # ------------------------------------------------------------------
    def onAnimateEndEvent(self, event):
        self._step += 1
        # If we exceed the pre-allocated buffer, grow it.  Should not happen
        # under normal runs but is cheap insurance.
        if self._step >= self._step_idx.size:
            self._grow_buffers(extra=1024)

        # ---- Per-step lightweight metrics --------------------------------
        wall_dt = (time.perf_counter() - self._t_begin) if self._t_begin else np.nan

        try:
            sim_time = float(self.root_node.time.value)
        except Exception:
            sim_time = np.nan

        n_iters  = self._read_solver_field(self._iter_field)
        residual = self._read_solver_field(self._error_field)

        max_pen, n_active = self._read_max_pen_and_count()

        i = self._step - 1  # 0-based store index
        self._step_idx[i]     = self._step
        self._sim_time[i]     = sim_time
        self._t2_pos_m[i]     = self.actuator._t2_pos_m
        self._t2_angle_rad[i] = self.actuator._t2_angle_rad
        self._n_iters[i]      = n_iters
        self._residual[i]     = residual
        self._max_pen[i]      = max_pen
        self._n_active[i]     = n_active
        self._wall_dt_s[i]    = wall_dt

        # ---- Snapshot (heavy data) at scheduled steps --------------------
        if self.schedule.is_snapshot_step(self._step):
            ph = self.schedule.phase_of(self._step)
            label = f"{ph['phase']}/{ph['sub']}/incr={ph['incr']}"
            s_B, gap_n = self._read_snapshot_profile()
            self._snap_step.append(self._step)
            self._snap_phase.append(label)
            self._snap_t2_pos.append(self.actuator._t2_pos_m)
            self._snap_t2_ang.append(self.actuator._t2_angle_rad)
            self._snap_s_B.append(s_B)
            self._snap_gap_n.append(gap_n)
            print(f"[ComparisonRecorder] snapshot {len(self._snap_step)}/"
                  f"{self.schedule.total_snapshots} at step={self._step} "
                  f"({label}): n_active={n_active}, "
                  f"max_pen={max_pen:+.3e} m", flush=True)

        # ---- End-of-run write --------------------------------------------
        if self._step >= self.schedule.total_steps:
            self._write_npz()

    # ==================================================================
    #  Helpers
    # ==================================================================

    @staticmethod
    def _probe_field(obj, candidates):
        """Return the first attribute name in `candidates` that resolves to a
        Data field with a readable .value, or None."""
        for name in candidates:
            try:
                f = getattr(obj, name)
                _ = f.value     # touch to ensure readable
                return name
            except Exception:
                continue
        return None

    @staticmethod
    def _list_available_data(obj):
        """List Data-field-like attributes on a SOFA object for debugging."""
        out = []
        for name in dir(obj):
            if name.startswith('_'):
                continue
            try:
                f = getattr(obj, name)
                if hasattr(f, 'value'):
                    out.append(name)
            except Exception:
                pass
        return ', '.join(sorted(out))

    def _read_solver_field(self, name):
        if name is None:
            return np.nan
        try:
            return float(getattr(self.solver, name).value)
        except Exception:
            return np.nan

    # ------------------------------------------------------------------
    def _read_max_pen_and_count(self):
        """
        Read BCM `distances`, filter out kInvalidGap sentinels, and return:
          max_pen  -- max penetration depth (= max of -delta_n over active
                      pairs where delta_n < 0), or 0.0 if all active pairs
                      have delta_n >= 0; NaN if no readable distances.
          n_active -- number of active pairs.
        """
        try:
            dists = np.asarray(self.bcm.distances.value)
        except Exception:
            return np.nan, -1
        if dists.size == 0:
            return 0.0, 0
        delta_n = dists[:, 0]
        valid = np.flatnonzero(delta_n < INVALID_GAP_THRESHOLD)
        if valid.size == 0:
            return 0.0, 0
        active = delta_n[valid]
        max_pen = float(max(0.0, -active.min())) if active.min() < 0 else 0.0
        return max_pen, int(valid.size)

    # ------------------------------------------------------------------
    def _read_snapshot_profile(self):
        """
        Heavy snapshot: filter active pairs from BCM.distances, read Pc_B
        from contactMO (contactPoints mode -> index 2k+1), project Pc_B
        onto Tube_3 centerline polyline to get curvilinear abscissa s_B,
        sort by s_B for clean plotting.
        """
        try:
            dists = np.asarray(self.bcm.distances.value)
            pos = np.asarray(self.contact_mo.position.value)
            centers = np.asarray(self.t2_MO.position.value)[:, :3]
        except Exception:
            return np.empty(0), np.empty(0)

        if dists.size == 0:
            return np.empty(0), np.empty(0)

        delta_n = dists[:, 0]
        valid_k = np.flatnonzero(delta_n < INVALID_GAP_THRESHOLD)
        if valid_k.size == 0:
            return np.empty(0), np.empty(0)

        idx_B = 2 * valid_k + 1
        if idx_B.max() >= pos.shape[0]:
            return np.empty(0), np.empty(0)
        if centers.shape[0] != self.t2_curv_abs.shape[0]:
            return np.empty(0), np.empty(0)

        pcB = pos[idx_B]
        gaps = delta_n[valid_k]
        absc = np.array([
            self._project_to_polyline(p, centers, self.t2_curv_abs)
            for p in pcB
        ])
        order = np.argsort(absc)
        return absc[order].astype(np.float64), gaps[order].astype(np.float64)

    @staticmethod
    def _project_to_polyline(point, centers, abscissae):
        """Lifted from LiveContactMonitor._project_to_polyline."""
        best_d2 = np.inf
        best_abs = float(abscissae[0])
        N = len(centers)
        for i in range(N - 1):
            C0 = centers[i]
            C1 = centers[i + 1]
            seg = C1 - C0
            seg_l2 = float(seg @ seg)
            if seg_l2 < 1e-20:
                continue
            t = float((point - C0) @ seg) / seg_l2
            t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
            proj = C0 + t * seg
            d2 = float((point - proj) @ (point - proj))
            if d2 < best_d2:
                best_d2 = d2
                best_abs = float(abscissae[i] + t * (abscissae[i + 1] - abscissae[i]))
        return best_abs

    # ------------------------------------------------------------------
    def _grow_buffers(self, extra):
        def _grow(arr, fill):
            new = np.full(arr.size + extra, fill, dtype=arr.dtype)
            new[:arr.size] = arr
            return new
        self._step_idx     = _grow(self._step_idx, -1)
        self._sim_time     = _grow(self._sim_time, np.nan)
        self._t2_pos_m     = _grow(self._t2_pos_m, np.nan)
        self._t2_angle_rad = _grow(self._t2_angle_rad, np.nan)
        self._n_iters      = _grow(self._n_iters, np.nan)
        self._residual     = _grow(self._residual, np.nan)
        self._max_pen      = _grow(self._max_pen, np.nan)
        self._n_active     = _grow(self._n_active, -1)
        self._wall_dt_s    = _grow(self._wall_dt_s, np.nan)

    # ------------------------------------------------------------------
    def _write_npz(self):
        if self._written:
            return
        self._written = True
        T = self._step  # actual number of steps recorded

        # Truncate per-step buffers to actual length.
        sl = slice(0, T)

        # Pack variable-length per-snapshot arrays as object arrays.
        snap_s_B   = np.empty(len(self._snap_s_B),  dtype=object)
        snap_gap_n = np.empty(len(self._snap_gap_n), dtype=object)
        for i, (s, g) in enumerate(zip(self._snap_s_B, self._snap_gap_n)):
            snap_s_B[i]   = np.asarray(s,  dtype=np.float64)
            snap_gap_n[i] = np.asarray(g, dtype=np.float64)

        os.makedirs(os.path.dirname(self.out_path) or '.', exist_ok=True)
        np.savez_compressed(
            self.out_path,
            mode_label=np.array(self.mode_label),
            # Schedule (so the post-processor can reconstruct phases)
            sched_init_steps          = np.int64(self.schedule.init_steps),
            sched_hold_steps          = np.int64(self.schedule.hold_steps),
            sched_trans_n_increments  = np.int64(self.schedule.trans_n_increments),
            sched_trans_steps_per_incr= np.int64(self.schedule.trans_steps_per_incr),
            sched_rot_n_increments    = np.int64(self.schedule.rot_n_increments),
            sched_rot_steps_per_incr  = np.int64(self.schedule.rot_steps_per_incr),
            # Per-step
            step_idx     = self._step_idx[sl],
            sim_time     = self._sim_time[sl],
            t2_pos_m     = self._t2_pos_m[sl],
            t2_angle_rad = self._t2_angle_rad[sl],
            n_iters      = self._n_iters[sl],
            residual     = self._residual[sl],
            max_pen      = self._max_pen[sl],
            n_active     = self._n_active[sl],
            wall_dt_s    = self._wall_dt_s[sl],
            # Per-snapshot
            snap_step_idx   = np.asarray(self._snap_step,    dtype=np.int64),
            snap_phase      = np.asarray(self._snap_phase,   dtype=object),
            snap_t2_pos_m   = np.asarray(self._snap_t2_pos,  dtype=np.float64),
            snap_t2_angle_rad=np.asarray(self._snap_t2_ang,  dtype=np.float64),
            snap_s_B        = snap_s_B,
            snap_gap_n      = snap_gap_n,
        )
        print(f"[ComparisonRecorder] wrote {T} step rows + "
              f"{len(self._snap_step)} snapshots to {self.out_path}",
              flush=True)
