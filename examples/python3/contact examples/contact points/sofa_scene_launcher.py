"""
SOFA Scene Launcher
-------------------
Small dropdown GUI to open SOFA scenes in runSofa without auto-starting them.
The window remains open after each launch so multiple scenes can run in parallel.
Run with:   python sofa_launcher.py
"""

import tkinter as tk
from tkinter import ttk, messagebox
import subprocess
import os


# =============================================================================
# CONFIGURATION  ---  EDIT THIS SECTION
# =============================================================================

# ---- Path to the runSofa executable ----------------------------------------
# Windows example: r"C:\Users\zawada\SOFA\build\install\bin\runSofa.exe"
# Linux   example:  "/home/user/sofa/build/install/bin/runSofa"
RUNSOFA_PATH = r"C:\Users\zawada\SOFA\build\bin\Release\runSofa.exe"


# ---- Scenes shown in the dropdown ------------------------------------------
# Each entry is a tuple: (display title, absolute path to the .py or .scn file)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>            FILL IN YOUR 4 TITLES AND 4 PATHS BELOW                    <
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
SCENES = [
    # ( "<<< TITLE shown in the dropdown >>>" ,  r"<<< PATH to scene file >>>" ),

    ( "CPUC: Cantilivear Beam",  r"C:\Users\zawada\SOFA\plugins\Cosserat\examples\python3\contact examples\contact points\cantiliver_beam.py" ),
    ( "CPUC: Sliding Cantilivear Beam",  r"C:\Users\zawada\SOFA\plugins\Cosserat\examples\python3\contact examples\contact points\sliding_cantilivear_beam.py" ),
    ( "CPUC: Straight Concentric Beam",  r"C:\Users\zawada\SOFA\plugins\Cosserat\examples\python3\contact examples\contact points\concentric_straight_beams.py" ),
    ( "CPUC: Two tube CTR",  r"C:\Users\zawada\SOFA\plugins\Cosserat\examples\python3\contact examples\contact points\ctr\ctr_two_tubes.py" ),
]


# ---- Optional extra args passed to runSofa ---------------------------------
# Leave empty for "loaded but paused" behaviour. Do NOT add "-a" / "--start"
# unless you want the simulation to auto-play on open.
EXTRA_RUNSOFA_ARGS = ["-l", "SofaPython3"]       # e.g. ["-g", "qglviewer"]


# =============================================================================
# Launcher GUI
# =============================================================================

class SofaSceneLauncher(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SOFA Scene Launcher")
        self.geometry("560x150")
        self.resizable(False, False)

        # Keep handles to spawned processes so they're not GC'd
        # and so we can report how many are still running.
        self._processes = []

        # ---- header ----
        ttk.Label(
            self,
            text="Select a scene to open in runSofa (loads paused):",
            font=("Segoe UI", 10),
        ).pack(pady=(12, 6))

        # ---- dropdown + button ----
        row = ttk.Frame(self)
        row.pack(padx=12, pady=4, fill="x")

        titles = [t for t, _ in SCENES]
        self.combo = ttk.Combobox(row, values=titles, state="readonly", width=50)
        if titles:
            self.combo.current(0)
        self.combo.pack(side="left", fill="x", expand=True, padx=(0, 8))
        # Double-click an item in the open dropdown to launch immediately:
        self.combo.bind("<<ComboboxSelected>>", lambda _e: None)

        ttk.Button(row, text="Open", command=self.open_selected).pack(side="left")

        # ---- status bar ----
        self.status = ttk.Label(self, text="Ready.", foreground="#444", anchor="w")
        self.status.pack(fill="x", padx=12, pady=(10, 12))

    # -------------------------------------------------------------------------
    def open_selected(self):
        idx = self.combo.current()
        if idx < 0:
            messagebox.showwarning("No scene selected",
                                   "Please pick a scene from the dropdown first.")
            return

        title, scene_path = SCENES[idx]

        # --- sanity checks ---
        if not os.path.isfile(RUNSOFA_PATH):
            messagebox.showerror(
                "runSofa not found",
                f"Could not find runSofa at:\n  {RUNSOFA_PATH}\n\n"
                "Edit RUNSOFA_PATH at the top of this script.")
            return

        if not os.path.isfile(scene_path):
            messagebox.showerror(
                "Scene file not found",
                f"Could not find scene file:\n  {scene_path}\n\n"
                "Edit the SCENES list at the top of this script.")
            return

        # --- launch (non-blocking) ---
        cmd = [RUNSOFA_PATH, scene_path, *EXTRA_RUNSOFA_ARGS]
        try:
            proc = subprocess.Popen(cmd)
            self._processes.append(proc)
            running = sum(1 for p in self._processes if p.poll() is None)
            self.status.config(
                text=f"Opened: {title}   (PID {proc.pid})   |   "
                     f"runSofa instances running: {running}")
        except Exception as exc:
            messagebox.showerror("Failed to launch runSofa", str(exc))


if __name__ == "__main__":
    SofaSceneLauncher().mainloop()