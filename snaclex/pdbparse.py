"""Minimal, dependency-free PDB-format parser.

Reads ATOM/HETATM records from the first model, classifies them into protein
atoms, water, monatomic ions/elements, and organic ligands, and groups hetero
components so they can be analyzed and selected in the UI.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "SOL"}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",  # selenomethionine + rare encoded residues
}

# Standard nucleotide residue names (DNA: DA/DC/DG/DT; RNA: A/C/G/U) plus a few
# common variants. Nucleotides arrive as ATOM records but are not amino acids, so
# without this set they were silently dropped from analysis (kept only for display).
NUCLEOTIDE_RESIDUES = {
    "DA", "DC", "DG", "DT", "DU", "DI",          # deoxyribonucleotides
    "A", "C", "G", "U", "I",                       # ribonucleotides
    "N",                                            # unknown nucleotide
    "5MC", "5MU", "1MA", "7MG", "2MG", "M2G", "OMC", "OMG", "PSU", "H2U",
}

# Common biologically relevant metals / monatomic ions seen as HETATM.
METAL_ELEMENTS = {
    "NA", "K", "MG", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "MO", "W",
    "CD", "HG", "PT", "AU", "AG", "PB", "BA", "SR", "CS", "RB", "LI", "AL",
    "V", "CR", "SN",
}
HALIDE_ELEMENTS = {"CL", "BR", "F", "I"}


@dataclass
class Atom:
    serial: int
    name: str
    res_name: str
    chain: str
    res_seq: int
    icode: str
    x: float
    y: float
    z: float
    element: str
    is_hetero: bool
    occupancy: float = 1.0
    bfactor: float = 0.0

    @property
    def key(self) -> tuple:
        return (self.chain, self.res_seq, self.icode)


@dataclass
class Component:
    """A grouped hetero component (ligand, ion, cofactor, etc.)."""

    res_name: str
    chain: str
    res_seq: int
    icode: str
    atoms: list[Atom] = field(default_factory=list)

    @property
    def label(self) -> str:
        seq = f"{self.res_seq}{self.icode}".strip()
        return f"{self.res_name} {self.chain}{seq}".strip()

    @property
    def kind(self) -> str:
        heavy = [a for a in self.atoms if a.element != "H"]
        if len(heavy) == 1:
            el = heavy[0].element
            if el in METAL_ELEMENTS:
                return "metal"
            if el in HALIDE_ELEMENTS:
                return "ion"
            return "ion"
        return "ligand"


@dataclass
class Structure:
    atoms: list[Atom]
    protein_atoms: list[Atom]
    components: list[Component]
    chains: list[str]
    nucleic_atoms: list[Atom] = field(default_factory=list)
    nucleic_chains: list[str] = field(default_factory=list)

    @property
    def ligand_components(self) -> list[Component]:
        return [c for c in self.components if c.kind == "ligand"]


def _f(s: str, default: float = 0.0) -> float:
    try:
        return float(s)
    except ValueError:
        return default


def parse_pdb(text: str) -> Structure:
    atoms: list[Atom] = []
    seen_alt: dict[tuple, str] = {}
    in_first_model = True

    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            in_first_model = False
            continue
        if rec == "MODEL":
            # Only keep the first model.
            if atoms:
                in_first_model = False
            continue
        if not in_first_model:
            continue
        if rec not in ("ATOM", "HETATM"):
            continue
        if len(line) < 54:
            continue

        alt_loc = line[16].strip()
        name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21].strip() or "A"
        res_seq = int(_f(line[22:26], 0))
        icode = line[26].strip()

        # Keep only one alternate location per atom (the first encountered).
        alt_key = (chain, res_seq, icode, name)
        if alt_loc:
            if alt_key in seen_alt and seen_alt[alt_key] != alt_loc:
                continue
            seen_alt[alt_key] = alt_loc

        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            # Fall back: derive from atom name.
            element = "".join(c for c in name if c.isalpha())[:2].upper()

        atoms.append(
            Atom(
                serial=int(_f(line[6:11], 0)),
                name=name,
                res_name=res_name,
                chain=chain,
                res_seq=res_seq,
                icode=icode,
                x=_f(line[30:38]),
                y=_f(line[38:46]),
                z=_f(line[46:54]),
                element=element,
                is_hetero=(rec == "HETATM"),
                occupancy=_f(line[54:60], 1.0) if len(line) >= 60 else 1.0,
                bfactor=_f(line[60:66], 0.0) if len(line) >= 66 else 0.0,
            )
        )

    return _assemble(atoms)


def subset_chain(structure: Structure, chain: str) -> Structure:
    """Return a new Structure containing only the atoms of one chain.

    Lets large multi-chain assemblies (which exceed PDB-format / interactive
    limits) be analyzed one chain at a time.
    """
    return _assemble([a for a in structure.atoms if a.chain == chain])


def _assemble(atoms: list[Atom]) -> Structure:
    """Group a flat atom list into protein atoms, nucleic atoms, hetero
    components, and chains."""
    protein_atoms = [a for a in atoms if not a.is_hetero and a.res_name in STANDARD_AA]
    nucleic_atoms = [
        a for a in atoms if not a.is_hetero and a.res_name in NUCLEOTIDE_RESIDUES
    ]

    comp_map: dict[tuple, Component] = {}
    for a in atoms:
        if not a.is_hetero:
            continue
        if a.res_name in WATER_NAMES:
            continue
        key = (a.res_name, a.chain, a.res_seq, a.icode)
        comp = comp_map.get(key)
        if comp is None:
            comp = Component(a.res_name, a.chain, a.res_seq, a.icode)
            comp_map[key] = comp
        comp.atoms.append(a)

    chains = sorted({a.chain for a in protein_atoms})
    nucleic_chains = sorted({a.chain for a in nucleic_atoms})
    components = sorted(
        comp_map.values(), key=lambda c: (c.chain, c.res_seq, c.res_name)
    )
    return Structure(
        atoms=atoms,
        protein_atoms=protein_atoms,
        components=components,
        chains=chains,
        nucleic_atoms=nucleic_atoms,
        nucleic_chains=nucleic_chains,
    )


def parse_structure(text: str) -> Structure:
    """Parse a structure file, auto-detecting PDB vs mmCIF/PDBx format.

    RCSB serves only mmCIF for entries with no legacy PDB-format file (large or
    newer structures), so the fetch layer falls back to ``.cif`` and we detect
    it here.
    """
    if "_atom_site." in text:
        return parse_mmcif(text)
    return parse_pdb(text)


_CIF_TOKEN = re.compile(r"'[^']*'|\"[^\"]*\"|\S+")


def _cif_tokens(line: str) -> list[str]:
    out = []
    for tok in _CIF_TOKEN.findall(line):
        if len(tok) >= 2 and tok[0] in "'\"" and tok[-1] == tok[0]:
            tok = tok[1:-1]
        out.append(tok)
    return out


def to_pdb(structure: Structure) -> str:
    """Serialize a parsed Structure back to PDB-format ATOM/HETATM text.

    Lets the 3Dmol viewer (which reads "pdb") render structures that arrived as
    mmCIF. PDB format caps serials at 99,999 and chain ids at one character, so
    this is only valid for structures within those limits (the size guard in the
    server enforces that before calling here).
    """
    lines = []
    for i, a in enumerate(structure.atoms, start=1):
        if i > 99999:
            break
        rec = "HETATM" if a.is_hetero else "ATOM  "
        name = (a.name or "")[:4]
        nm = name if len(name) >= 4 else " " + name
        chain = (a.chain or "A")[:1]
        res = (a.res_name or "")[:3]
        icode = ((a.icode or " ")[:1]) or " "
        el = (a.element or "")[:2].upper()
        lines.append(
            f"{rec}{i:>5} {nm:<4} {res:>3} {chain}{a.res_seq:>4}{icode}   "
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occupancy:6.2f}{a.bfactor:6.2f}"
            f"          {el:>2}"
        )
    lines.append("END")
    return "\n".join(lines)


def _cif_value(row: dict, *names, default: str = "") -> str:
    """First present, non-placeholder value among `names` ('.'/'?' mean unset)."""
    for n in names:
        v = row.get(n)
        if v not in (None, "", ".", "?"):
            return v
    return default


def parse_mmcif(text: str) -> Structure:
    """Parse the ``_atom_site`` loop of an mmCIF/PDBx file into a Structure.

    Handles the first model only and keeps a single alternate location per atom,
    mirroring ``parse_pdb``. Prefers auth_* identifiers (what users and the PDB
    format use) over label_* ones.
    """
    lines = text.splitlines()
    atoms: list[Atom] = []
    seen_alt: dict[tuple, str] = {}
    first_model: str | None = None

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.strip() != "loop_":
            i += 1
            continue
        # Collect the column tags that follow `loop_`.
        i += 1
        tags: list[str] = []
        while i < n and lines[i].lstrip().startswith("_"):
            tags.append(lines[i].strip())
            i += 1
        if not any(t.startswith("_atom_site.") for t in tags):
            continue  # some other loop; skip its tag block, scan onward
        col = {t: idx for idx, t in enumerate(tags)}

        def g(row, name):
            idx = col.get("_atom_site." + name)
            return row[idx] if idx is not None and idx < len(row) else None

        while i < n:
            raw = lines[i]
            s = raw.strip()
            if s == "" or s.startswith("#") or s.startswith("loop_") or s.startswith("_") or s.startswith("data_"):
                break
            i += 1
            row = _cif_tokens(raw)
            if len(row) < len(tags):
                continue

            model = g(row, "pdbx_PDB_model_num")
            if first_model is None:
                first_model = model
            elif model != first_model:
                continue

            d = {
                "group": g(row, "group_PDB") or "ATOM",
                "element": (g(row, "type_symbol") or "").strip().upper(),
                "name": _cif_value({"a": g(row, "auth_atom_id"), "l": g(row, "label_atom_id")}, "a", "l"),
                "alt": g(row, "label_alt_id") or "",
                "res_name": _cif_value({"a": g(row, "auth_comp_id"), "l": g(row, "label_comp_id")}, "a", "l"),
                "chain": _cif_value({"a": g(row, "auth_asym_id"), "l": g(row, "label_asym_id")}, "a", "l", default="A"),
                "seq": _cif_value({"a": g(row, "auth_seq_id"), "l": g(row, "label_seq_id")}, "a", "l", default="0"),
                "icode": g(row, "pdbx_PDB_ins_code") or "",
                "x": g(row, "Cartn_x"), "y": g(row, "Cartn_y"), "z": g(row, "Cartn_z"),
                "occ": g(row, "occupancy"), "b": g(row, "B_iso_or_equiv"),
            }
            if d["alt"] in (".", "?"):
                d["alt"] = ""
            if d["icode"] in (".", "?"):
                d["icode"] = ""
            name = d["name"].strip("'\"")
            res_seq = int(_f(d["seq"], 0))
            alt_key = (d["chain"], res_seq, d["icode"], name)
            if d["alt"]:
                if alt_key in seen_alt and seen_alt[alt_key] != d["alt"]:
                    continue
                seen_alt[alt_key] = d["alt"]

            element = d["element"] or "".join(c for c in name if c.isalpha())[:2].upper()
            atoms.append(
                Atom(
                    serial=int(_f(g(row, "id"), 0)),
                    name=name,
                    res_name=d["res_name"],
                    chain=d["chain"],
                    res_seq=res_seq,
                    icode=d["icode"],
                    x=_f(d["x"]), y=_f(d["y"]), z=_f(d["z"]),
                    element=element,
                    is_hetero=(d["group"] == "HETATM"),
                    occupancy=_f(d["occ"], 1.0),
                    bfactor=_f(d["b"], 0.0),
                )
            )
        break  # the atom_site loop is parsed; done

    return _assemble(atoms)
