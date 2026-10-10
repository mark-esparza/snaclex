"""Lightweight rigid-body molecular docking (pure Python, no external engine).

This is an AutoDock-style grid-map docker, simplified to run dependency-free:

  1. A scoring grid is precomputed over a box around the target pocket, with
     separate channels for steric packing, hydrogen bonding, and hydrophobic
     contact (mirroring AutoDock affinity maps).
  2. A real 3D ligand conformer (from PubChem) is searched as a rigid body via
     Monte-Carlo translation/rotation with simulated-annealing-style acceptance.
  3. The best pose's atoms are returned so its atomic interactions can be
     profiled with the same engine used for crystallographic ligands.

It is an approximate, research-only docker: the ligand is treated as rigid
(PubChem's single conformer), scoring is empirical and NOT calibrated to
kcal/mol, and there is no explicit solvent or full force field. Use it to
generate pose/interaction hypotheses, not affinities.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .pdbparse import Atom, Component, Structure

# Van der Waals radii (Angstrom).
_RVDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "B": 1.92, "SI": 2.10,
    "ZN": 1.39, "FE": 1.40, "MG": 1.73, "CA": 2.00, "NA": 2.27, "K": 2.75,
}
_DEFAULT_RVDW = 1.70

GRID_HALF = 10.0      # half-size of scoring box (A)
SPACING = 1.0         # grid spacing (A)
TRANS_HALF = 5.0      # how far the ligand center may wander from pocket center
STERIC_CUT = 6.5
HBOND_CUT = 3.6
HYDRO_CUT = 4.5
PROBE_R = 1.70        # generic carbon probe for the steric channel


def _rvdw(element: str) -> float:
    return _RVDW.get(element.upper(), _DEFAULT_RVDW)


@dataclass
class _RecAtom:
    x: float
    y: float
    z: float
    element: str
    is_acceptor: bool   # O -> acceptor
    is_donor: bool      # N -> donor
    is_carbon: bool


class _CellGrid:
    """Spatial hash for fast neighbor queries during grid construction."""

    def __init__(self, atoms, cell):
        self.cell = cell
        self.b = {}
        for a in atoms:
            k = (int(a.x // cell), int(a.y // cell), int(a.z // cell))
            self.b.setdefault(k, []).append(a)

    def near(self, x, y, z):
        cx, cy, cz = int(x // self.cell), int(y // self.cell), int(z // self.cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = self.b.get((cx + dx, cy + dy, cz + dz))
                    if bucket:
                        yield from bucket


@dataclass
class Grid:
    origin: tuple
    n: int
    spacing: float
    steric: list
    hb_don: list   # favorability for a ligand DONOR sitting here
    hb_acc: list   # favorability for a ligand ACCEPTOR sitting here
    hydro: list    # favorability for a ligand CARBON sitting here

    def _idx(self, i, j, k):
        return (i * self.n + j) * self.n + k

    def interp(self, x, y, z):
        """Trilinear interpolation of all channels; out-of-box => penalty."""
        ox, oy, oz = self.origin
        fx = (x - ox) / self.spacing
        fy = (y - oy) / self.spacing
        fz = (z - oz) / self.spacing
        i, j, k = int(fx), int(fy), int(fz)
        if i < 0 or j < 0 or k < 0 or i >= self.n - 1 or j >= self.n - 1 or k >= self.n - 1:
            return None  # signals out-of-box
        tx, ty, tz = fx - i, fy - j, fz - k
        s = d = a = h = 0.0
        for di in (0, 1):
            wx = tx if di else (1 - tx)
            for dj in (0, 1):
                wy = ty if dj else (1 - ty)
                for dk in (0, 1):
                    wz = tz if dk else (1 - tz)
                    w = wx * wy * wz
                    idx = self._idx(i + di, j + dj, k + dk)
                    s += w * self.steric[idx]
                    d += w * self.hb_don[idx]
                    a += w * self.hb_acc[idx]
                    h += w * self.hydro[idx]
        return s, d, a, h


def _build_receptor(structure: Structure, center, cutoff, extra_atoms=None) -> list:
    """Collect scoring atoms near `center`.

    `extra_atoms` adds non-protein atoms (cofactors, structural metals,
    catalytic waters) to the rigid receptor. By default only standard amino
    acids are scored, which silently removes FAD from a monoamine oxidase site,
    SAM from a methyltransferase, PLP from a decarboxylase, and so on -- so a
    cofactor-dependent site is scored as an empty cavity. Callers that know
    which heterocomponents are part of the site pass them here.
    """
    cx, cy, cz = center
    reach = GRID_HALF + cutoff
    rec = []
    for a in list(structure.protein_atoms) + list(extra_atoms or []):
        if a.element == "H":
            continue
        if abs(a.x - cx) > reach or abs(a.y - cy) > reach or abs(a.z - cz) > reach:
            continue
        el = a.element
        rec.append(
            _RecAtom(
                a.x, a.y, a.z, el,
                is_acceptor=(el == "O"),
                is_donor=(el == "N"),
                is_carbon=(el == "C"),
            )
        )
    return rec


def build_grid(structure: Structure, center, extra_atoms=None) -> Grid:
    """Precompute the scoring grid around `center`.

    Pass `extra_atoms` to include cofactors/metals as part of the rigid
    receptor; see :func:`_build_receptor`.
    """
    rec = _build_receptor(structure, center, STERIC_CUT, extra_atoms)
    cell = _CellGrid(rec, STERIC_CUT)
    n = int((2 * GRID_HALF) / SPACING) + 1
    ox = center[0] - GRID_HALF
    oy = center[1] - GRID_HALF
    oz = center[2] - GRID_HALF

    steric = [0.0] * (n * n * n)
    hb_don = [0.0] * (n * n * n)
    hb_acc = [0.0] * (n * n * n)
    hydro = [0.0] * (n * n * n)

    exp = math.exp
    sc2 = STERIC_CUT * STERIC_CUT
    hb2 = HBOND_CUT * HBOND_CUT
    hy2 = HYDRO_CUT * HYDRO_CUT

    for i in range(n):
        px = ox + i * SPACING
        for j in range(n):
            py = oy + j * SPACING
            base = (i * n + j) * n
            for k in range(n):
                pz = oz + k * SPACING
                s = d = a = h = 0.0
                for r in cell.near(px, py, pz):
                    ddx = px - r.x
                    ddy = py - r.y
                    ddz = pz - r.z
                    dist2 = ddx * ddx + ddy * ddy + ddz * ddz
                    if dist2 > sc2:
                        continue
                    dist = math.sqrt(dist2)
                    # steric (generic carbon probe)
                    delta = dist - (PROBE_R + _rvdw(r.element))
                    s += -0.8 * exp(-(delta / 0.5) ** 2)
                    s += -0.2 * exp(-((delta - 3.0) / 2.0) ** 2)
                    if delta < 0:
                        s += min(4.0 * delta * delta, 8.0)
                    # hydrogen bond channels
                    if dist2 <= hb2:
                        hd = dist - 2.9
                        bell = exp(-(hd / 0.5) ** 2)
                        if r.is_acceptor:
                            d += -1.2 * bell
                        if r.is_donor:
                            a += -1.2 * bell
                    # hydrophobic channel
                    if r.is_carbon and dist2 <= hy2:
                        h += -0.4 * exp(-((dist - 3.8) / 0.8) ** 2)
                steric[base + k] = s
                hb_don[base + k] = d
                hb_acc[base + k] = a
                hydro[base + k] = h

    return Grid((ox, oy, oz), n, SPACING, steric, hb_don, hb_acc, hydro)


# ---------------- ligand geometry ----------------

def _centroid(atoms):
    n = len(atoms)
    return (
        sum(a["x"] for a in atoms) / n,
        sum(a["y"] for a in atoms) / n,
        sum(a["z"] for a in atoms) / n,
    )


def _local_coords(atoms):
    cx, cy, cz = _centroid(atoms)
    return [(a["x"] - cx, a["y"] - cy, a["z"] - cz) for a in atoms]


def _rand_quat(rng):
    # Uniform random unit quaternion (Shoemake).
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    s1 = math.sqrt(1 - u1)
    s2 = math.sqrt(u1)
    return (
        s1 * math.sin(2 * math.pi * u2),
        s1 * math.cos(2 * math.pi * u2),
        s2 * math.sin(2 * math.pi * u3),
        s2 * math.cos(2 * math.pi * u3),
    )


def _quat_mul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    )


def _small_quat(rng, sigma):
    angle = rng.gauss(0, sigma)
    ax, ay, az = rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)
    norm = math.sqrt(ax * ax + ay * ay + az * az) or 1.0
    s = math.sin(angle / 2)
    return (math.cos(angle / 2), s * ax / norm, s * ay / norm, s * az / norm)


def _quat_matrix(q):
    w, x, y, z = q
    n = math.sqrt(w * w + x * x + y * y + z * z) or 1.0
    w, x, y, z = w / n, x / n, y / n, z / n
    return (
        1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
        2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
        2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y),
    )


def _place(local, q, center):
    m = _quat_matrix(q)
    cx, cy, cz = center
    out = []
    for (lx, ly, lz) in local:
        out.append(
            (
                m[0] * lx + m[1] * ly + m[2] * lz + cx,
                m[3] * lx + m[4] * ly + m[5] * lz + cy,
                m[6] * lx + m[7] * ly + m[8] * lz + cz,
            )
        )
    return out


def _score(coords, elements, grid: Grid) -> float:
    total = 0.0
    for (x, y, z), el in zip(coords, elements):
        vals = grid.interp(x, y, z)
        if vals is None:
            total += 6.0  # out-of-box penalty keeps ligand inside
            continue
        s, d, a, h = vals
        total += s
        if el == "N":
            total += d
        elif el == "O":
            total += a
        elif el == "C":
            total += h
    return total


# ---------------- search ----------------

def dock(structure: Structure, ligand_atoms: list, center,
         seeds: int = 220, mc_steps: int = 40, seed: int = 0,
         extra_atoms=None) -> dict:
    """Dock a rigid ligand into the pocket centered at `center` (builds grid)."""
    grid = build_grid(structure, center, extra_atoms)
    return dock_with_grid(grid, ligand_atoms, center, seeds, mc_steps, seed)


def dock_with_grid(grid: Grid, ligand_atoms: list, center,
                   seeds: int = 220, mc_steps: int = 40, seed: int = 0) -> dict:
    """Dock a rigid ligand against a prebuilt grid (reusable across ligands)."""
    if not ligand_atoms:
        raise ValueError("Ligand has no heavy atoms to dock")

    local = _local_coords(ligand_atoms)
    elements = [a["element"] for a in ligand_atoms]
    rng = random.Random(seed)
    cx, cy, cz = center

    best_score = float("inf")
    best_coords = None

    temperature = 2.0
    for _ in range(seeds):
        c = (
            cx + rng.uniform(-TRANS_HALF, TRANS_HALF),
            cy + rng.uniform(-TRANS_HALF, TRANS_HALF),
            cz + rng.uniform(-TRANS_HALF, TRANS_HALF),
        )
        q = _rand_quat(rng)
        coords = _place(local, q, c)
        cur = _score(coords, elements, grid)

        for step in range(mc_steps):
            sigma = 0.45 * (1 - step / mc_steps) + 0.1
            nc = (
                min(max(c[0] + rng.gauss(0, 0.6), cx - TRANS_HALF), cx + TRANS_HALF),
                min(max(c[1] + rng.gauss(0, 0.6), cy - TRANS_HALF), cy + TRANS_HALF),
                min(max(c[2] + rng.gauss(0, 0.6), cz - TRANS_HALF), cz + TRANS_HALF),
            )
            nq = _quat_mul(_small_quat(rng, sigma), q)
            ncoords = _place(local, nq, nc)
            ns = _score(ncoords, elements, grid)
            if ns < cur or rng.random() < math.exp(-(ns - cur) / temperature):
                c, q, coords, cur = nc, nq, ncoords, ns
                if cur < best_score:
                    best_score = cur
                    best_coords = coords

    return {
        "score": round(best_score, 2),
        "ligand_efficiency": round(best_score / len(elements), 3),
        "pose_coords": best_coords,
        "elements": elements,
        "center": [round(c, 2) for c in center],
        "box_half": GRID_HALF,
        "grid_spacing": SPACING,
        "translation_half": TRANS_HALF,
        "n_heavy_atoms": len(elements),
        "search": {"seeds": seeds, "mc_steps": mc_steps, "random_seed": seed},
    }


def pose_to_component(pose: dict, res_name: str = "LIG") -> Component:
    """Wrap a docked pose as a Component so the interaction profiler can run."""
    atoms = []
    el_counts: dict[str, int] = {}
    for (x, y, z), el in zip(pose["pose_coords"], pose["elements"]):
        el_counts[el] = el_counts.get(el, 0) + 1
        atoms.append(
            Atom(
                serial=len(atoms) + 1,
                name=f"{el}{el_counts[el]}",
                res_name=res_name,
                chain="X",
                res_seq=999,
                icode="",
                x=x, y=y, z=z,
                element=el,
                is_hetero=True,
            )
        )
    comp = Component(res_name, "X", 999, "", atoms)
    return comp


def pose_to_pdb(pose: dict, res_name: str = "LIG") -> str:
    """Serialize a docked pose to a minimal PDB block for the 3D viewer."""
    lines = []
    el_counts: dict[str, int] = {}
    for serial, ((x, y, z), el) in enumerate(
        zip(pose["pose_coords"], pose["elements"]), start=1
    ):
        el_counts[el] = el_counts.get(el, 0) + 1
        name = f"{el}{el_counts[el]}"[:4]
        lines.append(
            f"HETATM{serial:>5} {name:<4} {res_name:>3} X 999    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}"
        )
    lines.append("END")
    return "\n".join(lines)


def component_center(component: Component) -> tuple:
    n = len(component.atoms)
    return (
        sum(a.x for a in component.atoms) / n,
        sum(a.y for a in component.atoms) / n,
        sum(a.z for a in component.atoms) / n,
    )


def rmsd_to_reference(pose: dict, component: Component):
    """Nearest-atom RMSD (by element) between a docked pose and a reference.

    Used to validate *redocking*: when the docked chemical is the same molecule
    that is crystallographically bound, this approximates how close the predicted
    pose is to the experimental one. It is a nearest-neighbor RMSD (no exact atom
    correspondence), so it is a proxy, not a rigorous superposition RMSD.
    """
    ref = [(a.x, a.y, a.z, a.element) for a in component.atoms if a.element != "H"]
    if not ref:
        return None
    total = 0.0
    n = 0
    for (x, y, z), el in zip(pose["pose_coords"], pose["elements"]):
        best = None
        for (rx, ry, rz, rel) in ref:
            if rel != el:
                continue
            d2 = (x - rx) ** 2 + (y - ry) ** 2 + (z - rz) ** 2
            if best is None or d2 < best:
                best = d2
        if best is not None:
            total += best
            n += 1
    if n == 0:
        return None
    return round(math.sqrt(total / n), 2)
