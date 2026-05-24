"""Render py-SCORPIUS plots on the same data the R driver dumped."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))

import pyscorpius.plotting as pl
from ggplot2_py import ggsave


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)

    space = pd.read_csv(out_dir / "space.csv").to_numpy(dtype=np.float64)
    groups = pd.read_csv(out_dir / "groups.csv")["group"].tolist()
    path = pd.read_csv(out_dir / "path.csv").to_numpy(dtype=np.float64)
    time = pd.read_csv(out_dir / "time.csv")["time"].to_numpy(dtype=np.float64)

    p1 = pl.draw_trajectory_plot(space, progression_group=groups)
    ggsave(str(out_dir / "Py_traj_plain.png"), plot=p1, width=5, height=5, dpi=100)

    p2 = pl.draw_trajectory_plot(space, progression_group=groups, path=path)
    ggsave(str(out_dir / "Py_traj_path.png"), plot=p2, width=5, height=5, dpi=100)

    p3 = pl.draw_trajectory_plot(space, progression_group=time)
    ggsave(str(out_dir / "Py_traj_numeric.png"), plot=p3, width=5, height=5, dpi=100)

    print("done")


if __name__ == "__main__":
    main()
