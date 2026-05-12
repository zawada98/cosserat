#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_comparison.py
================================================================================
Driver: launches `runSofa --start -g batch compare_scene.py` once per mode
(cpulc, then feeder), sequentially, capturing stdout/stderr to log files.

Each child process runs in its own subprocess so they don't share Python
state.  Sequential (not parallel) so the wall-clock measurements aren't
contaminated by CPU contention between the two runs.

Usage:
    python run_comparison.py
    python run_comparison.py --modes cpulc            # just one mode
    python run_comparison.py --runsofa /path/to/runSofa
    python run_comparison.py --schedule init_only     # short dry pass
                                                      # (NOT YET IMPLEMENTED
                                                      #  in compare_scene.py;
                                                      #  reserved for future)

Outputs:
    ./comparison_cpulc.npz
    ./comparison_feeder.npz
    ./run_<mode>.log
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def find_runsofa(explicit: str | None) -> str:
    if explicit:
        return explicit
    found = shutil.which('runSofa')
    if found:
        return found
    raise SystemExit(
        "runSofa not found on PATH.  Pass --runsofa /full/path/to/runSofa."
    )


def run_one(mode: str, runsofa: str, scene: Path,
            out_dir: Path, source_dir: Path | None,
            n_iters: int) -> dict:
    out_path = out_dir / f"comparison_{mode}.npz"
    log_path = out_dir / f"run_{mode}.log"

    env = os.environ.copy()
    env['CTR_COMPARE_MODE'] = mode
    env['CTR_COMPARE_OUT'] = str(out_path.resolve())
    if source_dir is not None:
        env['CTR_COMPARE_SOURCE'] = str(source_dir.resolve())

    cmd = [runsofa, '--start', '-g', 'batch',
           '-l', 'SofaPython3',
           '-n', str(n_iters),
           str(scene.resolve())]

    print(f"\n=========================================================="
          f"\n  RUN  mode={mode}"
          f"\n  cmd: {' '.join(cmd)}"
          f"\n  log: {log_path}"
          f"\n  out: {out_path}"
          f"\n==========================================================\n",
          flush=True)

    t0 = time.perf_counter()
    with open(log_path, 'w') as logf:
        proc = subprocess.run(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - t0

    print(f"[run_comparison] mode={mode} returncode={proc.returncode} "
          f"wall={wall:.1f}s", flush=True)
    return {
        'mode': mode, 'rc': proc.returncode, 'wall_s': wall,
        'out': out_path, 'log': log_path,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--modes', default='cpulc,feeder',
                    help="Comma-separated modes to run (default both).")
    ap.add_argument('--runsofa', default=None,
                    help="Full path to runSofa binary.")
    ap.add_argument('--scene', default=None,
                    help="Path to compare_scene.py (default: alongside this script).")
    ap.add_argument('--out-dir', default='.',
                    help="Output directory for .npz and .log files.")
    ap.add_argument('--source-dir', default=None,
                    help="Directory containing ctr_two_tubes.py (passed via "
                         "CTR_COMPARE_SOURCE env var if set).")
    ap.add_argument('--n-iters', type=int, default=30_000,
                    help="runSofa batch iteration count (default 30000, "
                         "covers the full 173k-step schedule with margin).")
    args = ap.parse_args()

    runsofa = find_runsofa(args.runsofa)
    here = Path(__file__).resolve().parent
    scene = Path(args.scene).resolve() if args.scene else (here / 'compare_scene.py')
    if not scene.exists():
        raise SystemExit(f"Scene file not found: {scene}")

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    source_dir = Path(args.source_dir).resolve() if args.source_dir else None

    modes = [m.strip().lower() for m in args.modes.split(',') if m.strip()]
    for m in modes:
        if m not in ('cpulc', 'feeder'):
            raise SystemExit(f"Unknown mode {m!r}; expected cpulc or feeder.")

    results = []
    for m in modes:
        results.append(run_one(m, runsofa, scene, out_dir, source_dir, args.n_iters))

    print("\n=== Summary ===")
    for r in results:
        ok = 'OK' if r['rc'] == 0 else f"FAIL (rc={r['rc']})"
        print(f"  mode={r['mode']:6}  {ok:14}  wall={r['wall_s']:8.1f}s  "
              f"out={r['out'].name}")
    if any(r['rc'] != 0 for r in results):
        sys.exit(1)


if __name__ == '__main__':
    main()
