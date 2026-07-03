"""Large-scale motion via Normal Mode Analysis (Anisotropic Network Model).

This is analytic linear algebra on the Cα coordinates you already have — NOT
molecular dynamics, no force field, no GPU. Each residue is a bead; beads within
a cutoff are connected by springs; we build the 3N×3N Hessian and solve for the
lowest-frequency non-trivial modes (skipping the 6 trivial rigid-body modes).
Low-frequency modes are the large-scale domain/hinge "breathing" motions.

Compute budget (the app is hosted on a small tier): a full 3N×3N eigensolver in
pure Python is too slow for big proteins, so we **coarse-grain** to at most
``MAX_BEADS`` beads (documented cap) and extract only the lowest modes via
shift-invert subspace iteration — Cholesky-factorize (H + εI) once, then iterate
with its inverse so the low-frequency end converges quickly. Everything stays
within a synchronous request.
"""

from __future__ import annotations

import math
import random

# Coarse-graining / method parameters (documented cap).
MAX_BEADS = 60          # residues are grouped into at most this many beads
CUTOFF_A = 16.0         # spring cutoff between beads (Å); coarse beads -> larger
GAMMA = 1.0             # uniform spring constant (arbitrary units)
N_MODES_DEFAULT = 6     # non-trivial low-frequency modes reported
_N_TRIVIAL = 6          # rigid-body (near-zero) modes to skip
_SHIFT = 1e-6           # diagonal shift making (H+εI) positive-definite
_SUBSPACE_ITERS = 28    # shift-invert subspace-iteration sweeps


def _ca_atoms(structure):
    """Representative Cα (fallback: first atom) per residue, in chain/seq order."""
    seen = {}
    order = []
    for a in structure.protein_atoms:
        key = (a.chain, a.res_seq)
        if a.name == "CA":
            seen[key] = a
        elif key not in seen:
            seen[key] = a
        if key not in order:
            order.append(key)
    return [(k, seen[k]) for k in order if k in seen]


def _beads(structure, cap):
    """Coarse-grain residues into <= cap beads (centroids of contiguous groups).

    Returns (positions, members) where members[b] is the list of (chain, res_seq)
    residues represented by bead b — used to spread each bead's mode vector back
    onto all its residues for animation.
    """
    cas = _ca_atoms(structure)
    n = len(cas)
    if n == 0:
        return [], []
    if n <= cap:
        positions = [(a.x, a.y, a.z) for _k, a in cas]
        members = [[k] for k, _a in cas]
        return positions, members

    # Split into `cap` contiguous groups; each bead = centroid of its group.
    positions, members = [], []
    per = n / cap
    for b in range(cap):
        lo = int(round(b * per))
        hi = int(round((b + 1) * per))
        hi = max(hi, lo + 1)
        group = cas[lo:hi]
        if not group:
            continue
        cx = sum(a.x for _k, a in group) / len(group)
        cy = sum(a.y for _k, a in group) / len(group)
        cz = sum(a.z for _k, a in group) / len(group)
        positions.append((cx, cy, cz))
        members.append([k for k, _a in group])
    return positions, members


def _hessian(positions, cutoff, gamma):
    """Dense 3N×3N ANM Hessian as a flat row-major list."""
    n = len(positions)
    dim = 3 * n
    H = [0.0] * (dim * dim)
    c2 = cutoff * cutoff
    for i in range(n):
        xi, yi, zi = positions[i]
        for j in range(i + 1, n):
            xj, yj, zj = positions[j]
            dx, dy, dz = xj - xi, yj - yi, zj - zi
            d2 = dx * dx + dy * dy + dz * dz
            if d2 > c2 or d2 == 0.0:
                continue
            k = gamma / d2
            # 3×3 super-element blocks (i,j),(j,i) and diagonal accumulation.
            bxx, bxy, bxz = k * dx * dx, k * dx * dy, k * dx * dz
            byy, byz, bzz = k * dy * dy, k * dy * dz, k * dz * dz
            block = ((bxx, bxy, bxz), (bxy, byy, byz), (bxz, byz, bzz))
            ib, jb = 3 * i, 3 * j
            for p in range(3):
                for q in range(3):
                    v = block[p][q]
                    # off-diagonal i-j (negative)
                    H[(ib + p) * dim + (jb + q)] -= v
                    H[(jb + p) * dim + (ib + q)] -= v
                    # diagonal blocks (positive)
                    H[(ib + p) * dim + (ib + q)] += v
                    H[(jb + p) * dim + (jb + q)] += v
    return H, dim


def _cholesky(A, n):
    """In-place-ish Cholesky of SPD flat matrix A (n×n) -> lower L (flat)."""
    L = [0.0] * (n * n)
    for i in range(n):
        Li = i * n
        for j in range(i + 1):
            Lj = j * n
            s = A[Li + j]
            for k in range(j):
                s -= L[Li + k] * L[Lj + k]
            if i == j:
                if s <= 0:
                    s = 1e-12
                L[Li + j] = math.sqrt(s)
            else:
                L[Li + j] = s / L[Lj + j]
    return L


def _solve_chol(L, n, b):
    """Solve (L Lᵀ) x = b via forward+back substitution."""
    y = [0.0] * n
    for i in range(n):
        Li = i * n
        s = b[i]
        for k in range(i):
            s -= L[Li + k] * y[k]
        y[i] = s / L[Li + i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for k in range(i + 1, n):
            s -= L[k * n + i] * x[k]
        x[i] = s / L[i * n + i]
    return x


def _matvec(A, n, x):
    out = [0.0] * n
    for i in range(n):
        Ai = i * n
        s = 0.0
        for k in range(n):
            s += A[Ai + k] * x[k]
        out[i] = s
    return out


def _orthonormalize(cols, n):
    """Modified Gram-Schmidt on a list of column vectors (in place)."""
    for i in range(len(cols)):
        vi = cols[i]
        for j in range(i):
            vj = cols[j]
            dot = sum(vi[k] * vj[k] for k in range(n))
            for k in range(n):
                vi[k] -= dot * vj[k]
        norm = math.sqrt(sum(v * v for v in vi)) or 1.0
        cols[i] = [v / norm for v in vi]


def _jacobi(A, n, sweeps=60):
    """Symmetric eigensolver (Jacobi) for small dense flat matrices."""
    a = list(A)
    V = [1.0 if i == j else 0.0 for i in range(n) for j in range(n)]
    for _ in range(sweeps):
        off = 0.0
        for p in range(n):
            for q in range(p + 1, n):
                off += a[p * n + q] ** 2
        if off < 1e-18:
            break
        for p in range(n):
            for q in range(p + 1, n):
                apq = a[p * n + q]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = a[p * n + p], a[q * n + q]
                phi = 0.5 * math.atan2(2 * apq, aqq - app)
                c, s = math.cos(phi), math.sin(phi)
                for k in range(n):
                    akp, akq = a[k * n + p], a[k * n + q]
                    a[k * n + p] = c * akp - s * akq
                    a[k * n + q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p * n + k], a[q * n + k]
                    a[p * n + k] = c * apk - s * aqk
                    a[q * n + k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k * n + p], V[k * n + q]
                    V[k * n + p] = c * vkp - s * vkq
                    V[k * n + q] = s * vkp + c * vkq
    eig = [a[i * n + i] for i in range(n)]
    vecs = [[V[r * n + c] for r in range(n)] for c in range(n)]  # vecs[c] = column c
    return eig, vecs


def _rigid_basis(positions):
    """The 6 rigid-body modes (3 translations + 3 rotations about the centroid).

    These are the analytic null space of the ANM Hessian. Projecting them out of
    the subspace iteration confines the search to *internal* motions, so the
    lowest non-trivial modes converge fast (instead of the shift-invert step
    being swamped by the rigid modes it amplifies most).
    """
    n = len(positions)
    dim = 3 * n
    cx = sum(p[0] for p in positions) / n
    cy = sum(p[1] for p in positions) / n
    cz = sum(p[2] for p in positions) / n
    basis = [[0.0] * dim for _ in range(6)]
    for i, (x, y, z) in enumerate(positions):
        rx, ry, rz = x - cx, y - cy, z - cz
        basis[0][3 * i] = 1.0                       # translate x
        basis[1][3 * i + 1] = 1.0                   # translate y
        basis[2][3 * i + 2] = 1.0                   # translate z
        basis[3][3 * i + 1] = -rz; basis[3][3 * i + 2] = ry   # rotate about x
        basis[4][3 * i] = rz;      basis[4][3 * i + 2] = -rx  # rotate about y
        basis[5][3 * i] = -ry;     basis[5][3 * i + 1] = rx   # rotate about z
    _orthonormalize(basis, dim)
    return basis


def _project_out(vec, basis, dim):
    for b in basis:
        dot = 0.0
        for k in range(dim):
            dot += vec[k] * b[k]
        for k in range(dim):
            vec[k] -= dot * b[k]


def _lowest_modes(H, dim, positions, n_modes):
    """Lowest n_modes *non-trivial* eigenpairs of H (rigid modes projected out)."""
    A = list(H)
    for i in range(dim):
        A[i * dim + i] += _SHIFT
    L = _cholesky(A, dim)
    rigid = _rigid_basis(positions)

    p = min(dim - 6, n_modes + 3)
    rng = random.Random(0)
    X = [[rng.gauss(0, 1) for _ in range(dim)] for _ in range(p)]
    for col in X:
        _project_out(col, rigid, dim)
    _orthonormalize(X, dim)
    for _ in range(_SUBSPACE_ITERS):
        X = [_solve_chol(L, dim, col) for col in X]
        for col in X:
            _project_out(col, rigid, dim)   # stay orthogonal to the rigid space
        _orthonormalize(X, dim)

    # Rayleigh-Ritz on H within the converged (rigid-free) subspace.
    HX = [_matvec(H, dim, col) for col in X]
    B = [0.0] * (p * p)
    for i in range(p):
        for j in range(i, p):
            v = sum(X[i][k] * HX[j][k] for k in range(dim))
            B[i * p + j] = B[j * p + i] = v
    ritz_vals, ritz_vecs = _jacobi(B, p)
    order = sorted(range(p), key=lambda idx: ritz_vals[idx])

    modes = []
    for idx in order[:n_modes]:
        lam = ritz_vals[idx]
        coeff = ritz_vecs[idx]
        vec = [0.0] * dim
        for c in range(p):
            cc = coeff[c]
            Xc = X[c]
            for k in range(dim):
                vec[k] += cc * Xc[k]
        modes.append((max(lam, 0.0), vec))
    return modes


def _collectivity(res_disp):
    """Fraction-of-structure-that-moves metric (κ), 0..1."""
    mags2 = [d[0] ** 2 + d[1] ** 2 + d[2] ** 2 for d in res_disp]
    total = sum(mags2)
    if total <= 0:
        return 0.0
    ent = 0.0
    for m in mags2:
        p = m / total
        if p > 0:
            ent -= p * math.log(p)
    return round(math.exp(ent) / len(mags2), 3)


def analyze(structure, n_modes=N_MODES_DEFAULT, cap=MAX_BEADS, cutoff=CUTOFF_A):
    """Compute the lowest non-trivial ANM modes for a structure.

    Returns a dict with per-residue displacement vectors for each mode (spread
    from beads onto member residues), eigenvalues, relative frequencies, and
    collectivity. ``available`` is False for structures too small to analyze.
    """
    positions, members = _beads(structure, cap)
    n = len(positions)
    if n < 6:
        return {"available": False,
                "reason": "Too few residues for normal-mode analysis."}

    H, dim = _hessian(positions, cutoff, GAMMA)
    # Rigid-body modes are projected out inside the solver, so these are already
    # the lowest-frequency *internal* (non-trivial) modes.
    nontrivial = _lowest_modes(H, dim, positions, n_modes)

    # Compact payload: each bead carries its member residues once; each mode
    # carries one 3-vector per bead. The client expands bead -> residues for
    # animation and contact-strain (many residues share a bead after coarsening).
    bead_members = [[[chain, res_seq] for (chain, res_seq) in members[b]]
                    for b in range(n)]

    modes = []
    for mi, (lam, vec) in enumerate(nontrivial):
        bead_vecs = [(vec[3 * b], vec[3 * b + 1], vec[3 * b + 2]) for b in range(n)]
        maxmag = max(math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2) for v in bead_vecs) or 1.0
        bead_disp = [[round(v[0] / maxmag, 4), round(v[1] / maxmag, 4),
                      round(v[2] / maxmag, 4)] for v in bead_vecs]
        res_disp = [d for b in range(n) for d in [bead_disp[b]] * len(members[b])]
        modes.append({
            "index": mi,
            "eigenvalue": round(lam, 6),
            "frequency": round(math.sqrt(lam), 4),
            "collectivity": _collectivity(res_disp),
            "bead_disp": bead_disp,
        })

    # Relative frequency (1 = the highest of the reported low modes).
    fmax = max((m["frequency"] for m in modes), default=1.0) or 1.0
    for m in modes:
        m["relative_frequency"] = round(m["frequency"] / fmax, 3)

    residue_count = sum(len(g) for g in members)
    return {
        "available": True,
        "n_beads": n,
        "residue_count": residue_count,
        "coarse_grained": residue_count > n,
        "cutoff_A": cutoff,
        "max_beads": cap,
        "n_modes": len(modes),
        "bead_members": bead_members,
        "modes": modes,
    }
