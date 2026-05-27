# -*- coding: utf-8 -*-
"""Live 3D shape monitor for the protruded part of Tube_3."""

import math

import numpy as np
import Sofa.Core


def _q_conj(q):
    return np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)


def _q_mul(a, b):
    ax, ay, az, aw = a
    bx, by, bz, bw = b
    return np.array([
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ], dtype=float)


def _q_normalized(q):
    n = float(np.linalg.norm(q))
    if n < 1e-15:
        return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
    return np.asarray(q, dtype=float) / n


def _q_rotate(q, v):
    q = _q_normalized(q)
    p = np.array([v[0], v[1], v[2], 0.0], dtype=float)
    return _q_mul(_q_mul(q, p), _q_conj(q))[:3]


def _q_from_kappa_z(ds, kappa):
    angle = float(kappa) * float(ds)
    return np.array([0.0, 0.0, math.sin(0.5 * angle), math.cos(0.5 * angle)],
                    dtype=float)


class ProtrudedShapeMonitor(Sofa.Core.Controller):
    """
    Push the actual protruded inner-tube centerline and its natural shape.

    The actual curve is read directly from Tube_3 FramesMO positions.  The
    natural curve is integrated from the same emergence frame using Tube_3's
    rest Cosserat strain DOFs, so it represents the shape the exposed part
    would take from precurvature alone with the same attachment pose.
    """

    def __init__(self,
                 gui_bridge,
                 t1_base_mo,
                 t2_base_mo,
                 t2_frames_mo,
                 t2_coss_mo,
                 t2_frame_curv_abs,
                 t2_sec_curv_abs,
                 t2_length,
                 t2_x_offset,
                 every_n_steps=20,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gui = gui_bridge
        self.t1_base_mo = t1_base_mo
        self.t2_base_mo = t2_base_mo
        self.t2_frames_mo = t2_frames_mo
        self.t2_coss_mo = t2_coss_mo
        self.t2_frame_curv_abs = np.asarray(t2_frame_curv_abs, dtype=float)
        self.t2_sec_curv_abs = np.asarray(t2_sec_curv_abs, dtype=float)
        self.t2_length = float(t2_length)
        self.t2_x_offset = float(t2_x_offset)
        self.every_n_steps = max(1, int(every_n_steps))
        self.step = 0

    def _rigid_x(self, mo):
        return float(mo.position.value[0][0])

    def onAnimateEndEvent(self, event):
        self.step += 1
        if self.step % self.every_n_steps:
            return

        protrusion = self._protrusion_length()
        if protrusion <= 1e-12:
            self._push_empty(protrusion)
            return

        frames = np.asarray(self.t2_frames_mo.position.value, dtype=float)
        if frames.ndim != 2 or frames.shape[1] < 7:
            self._push_empty(protrusion)
            return

        n_frames = min(frames.shape[0], self.t2_frame_curv_abs.shape[0])
        if n_frames == 0:
            self._push_empty(protrusion)
            return

        frames = frames[:n_frames]
        frame_s = self.t2_frame_curv_abs[:n_frames]
        start_s = max(0.0, self.t2_length - protrusion)
        ids = np.flatnonzero(frame_s >= start_s - 1e-12)
        if ids.size == 0:
            self._push_empty(protrusion)
            return

        anchor = int(ids[0])
        selected = ids
        actual_centers = frames[selected, :3]

        natural_centers = self._integrate_natural_from_anchor(
            frames[anchor, :3],
            frames[anchor, 3:7],
            frame_s[anchor],
            frame_s[selected])

        push = getattr(self.gui, 'push_protruded_shape_profile', None)
        if push is not None:
            push(self.step, actual_centers, natural_centers, protrusion)

    def _protrusion_length(self):
        t1_advance = self._rigid_x(self.t1_base_mo)
        t2_advance = self._rigid_x(self.t2_base_mo) - self.t2_x_offset
        return max(0.0, min(self.t2_length, t2_advance - t1_advance))

    def _push_empty(self, protrusion):
        push = getattr(self.gui, 'push_protruded_shape_profile', None)
        if push is not None:
            empty = np.empty((0, 3), dtype=float)
            push(self.step, empty, empty, protrusion)

    def _integrate_natural_from_anchor(self, anchor_pos, anchor_quat,
                                       anchor_s, targets_s):
        rests = np.asarray(self.t2_coss_mo.rest_position.value, dtype=float)
        n_sec = min(rests.shape[0], self.t2_sec_curv_abs.shape[0] - 1)
        if n_sec <= 0:
            return np.repeat(np.asarray(anchor_pos, dtype=float)[None, :],
                             len(targets_s), axis=0)

        sec_abs = self.t2_sec_curv_abs[:n_sec + 1]
        targets = np.asarray(targets_s, dtype=float)
        out = []
        pos = np.asarray(anchor_pos, dtype=float).copy()
        quat = _q_normalized(anchor_quat)
        current_s = float(anchor_s)

        for target_s in targets:
            target_s = float(max(current_s, target_s))
            while current_s < target_s - 1e-12:
                sec_i = int(np.searchsorted(sec_abs, current_s, side='right') - 1)
                sec_i = max(0, min(sec_i, n_sec - 1))
                sec_end = min(float(sec_abs[sec_i + 1]), target_s)
                ds = max(0.0, sec_end - current_s)
                kappa = float(rests[sec_i, 2])

                if abs(kappa) < 1e-12:
                    local_delta = np.array([ds, 0.0, 0.0], dtype=float)
                else:
                    local_delta = np.array([
                        math.sin(kappa * ds) / kappa,
                        (1.0 - math.cos(kappa * ds)) / kappa,
                        0.0,
                    ], dtype=float)

                pos += _q_rotate(quat, local_delta)
                quat = _q_normalized(_q_mul(quat, _q_from_kappa_z(ds, kappa)))
                current_s += ds

            out.append(pos.copy())

        return np.asarray(out, dtype=float)
