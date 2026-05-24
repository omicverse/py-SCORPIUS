"""Parity gate: pseudotime Pearson + low-dim space Procrustes vs R."""
import json, sys, subprocess
from pathlib import Path
import numpy as np
import pytest
import yaml

PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORT))
sys.path.insert(0, str(PORT.parent / "omicverse-rebuildr" / "engine"))
from parity_metrics import compute_parity, is_pass


@pytest.fixture(scope="session")
def manifest():
    return yaml.safe_load((PORT / "data" / "manifest.yaml").read_text())


@pytest.fixture(scope="session")
def parity_outputs():
    ref = PORT / "data" / "reference_output.json"
    cand = PORT / "data" / "candidate_output.json"
    if not ref.exists() or not cand.exists():
        pytest.skip("Reference / candidate outputs not generated — run "
                    "r_reference_driver.R and _run_candidate.py first.")
    return json.loads(ref.read_text()), json.loads(cand.read_text())


def test_pseudotime_parity(manifest, parity_outputs):
    r, p = parity_outputs
    r_pt = np.array(r["pseudotime"]); p_pt = np.array(p["pseudotime"])
    # Pseudotime direction is arbitrary — try both
    fwd = compute_parity(r_pt, p_pt, "ordinal")
    rev = compute_parity(r_pt, 1 - p_pt, "ordinal")
    best = max(fwd, rev)
    threshold = manifest["outputs"][0]["threshold"]
    assert best >= threshold, f"pseudotime Pearson {best:.4f} < {threshold}"


def test_low_dim_space_parity(manifest, parity_outputs):
    r, p = parity_outputs
    r_sp = np.array(r["space"]); p_sp = np.array(p["space"])
    proc = compute_parity(r_sp, p_sp, "embedding")
    threshold = manifest["outputs"][1]["threshold"]
    assert proc >= threshold, f"space Procrustes {proc:.4f} < {threshold}"
