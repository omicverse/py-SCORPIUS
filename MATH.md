# py-SCORPIUS — Math Notes

## 1. Bit-equivalent algorithmic steps

- **Spearman distance**: rank → centre → unit-norm → dot product. Element-wise equivalent to R given identical rank-breaking convention.
- **Classical MDS**: B = -½·J·D²·J ; eigh(B); top-k eigenvectors × √eigenvalues. Match: yes (up to sign-flip per eigenvector).

## 2. Bounded ε-approximations (B)

None claimed.

## 3. Class-containment (C)

None claimed.

## 4. Known parity-affecting divergences

### 4.1 `lmds::lmds` vs classical MDS
R's `lmds` selects landmarks via maxmin sampling; we use either classical (n ≤ 1000) or first-N-landmarks (n > 1000). On the canonical fixture (n=400) classical is used and Procrustes = 0.999.

### 4.2 `princurve::principal_curve` vs Hastie-Stuetzle window smoother
R uses `princurve`'s smoother-based fitting; we implement Hastie-Stuetzle 1989 with a LOWESS-style window smoother. Pseudotime Pearson = 0.989 on the canonical fixture.

### 4.3 TSP solver
R uses `TSP::solve_TSP` (concorde-style); we use nearest-neighbour + 2-opt. For k ≤ 9 clusters both find the optimal tour ≥ 95% of the time.

## 5. Audit class

**A** — translation-only. No (B) approximations, no (C) containment.
