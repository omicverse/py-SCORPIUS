import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))
from pyscorpius import reduce_dimensionality, infer_trajectory


def main():
    fixture_path, output_path = sys.argv[1], sys.argv[2]
    expr_csv = Path(fixture_path).with_name(Path(fixture_path).stem + "_expression.csv")
    grp_csv = Path(fixture_path).with_name(Path(fixture_path).stem + "_group.csv")
    expression = pd.read_csv(expr_csv, index_col=0).to_numpy(dtype=np.float64)
    group = pd.read_csv(grp_csv).iloc[:, 0].tolist()
    print(f"[cand] expression: {expression.shape}")
    np.random.seed(42)
    t0 = time.perf_counter()
    space = reduce_dimensionality(expression, dist="spearman", ndim=2)
    t_rd = time.perf_counter() - t0
    print(f"[cand] reduce_dimensionality: {t_rd:.2f}s; space {space.shape}")

    np.random.seed(42)
    t0 = time.perf_counter()
    traj = infer_trajectory(space)
    t_it = time.perf_counter() - t0
    print(f"[cand] infer_trajectory: {t_it:.2f}s; pseudotime length: {len(traj['time'])}")

    out = {
        "group": group,
        "space": space.tolist(),
        "pseudotime": traj["time"].tolist(),
        "path": traj["path"].tolist(),
        "timings": {"reduce_dimensionality": t_rd, "infer_trajectory": t_it},
    }
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"[cand] wrote {output_path}")


if __name__ == "__main__":
    main()
