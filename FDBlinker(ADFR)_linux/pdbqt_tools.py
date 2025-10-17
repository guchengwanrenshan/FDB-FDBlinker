#!/usr/bin/env python3
"""
pdbqt_tools.py — Utilities to merge/convert PDBQT and write PDB with CONECT,
including robust element typing and ARG hydrogen name normalization.

Changes in this revision:
- RENAMES non‑standard ARG guanidinium H names:
    1HH1→HH11, 2HH1→HH12, 1HH2→HH21, 2HH2→HH22
  (previous version incorrectly dropped these atoms in PDB output)
"""

from __future__ import annotations
import sys
from pathlib import Path
from math import sqrt
import os
import re
import shutil
import subprocess
from typing import List, Dict, Tuple, Iterable, Optional

# --- Basic periodic table (for element-name detection) ---
ELEMENTS = {
    # 1–36
    "H","HE","LI","BE","B","C","N","O","F","NE","NA","MG","AL","SI","P","S","CL","AR",
    "K","CA","SC","TI","V","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR",
    # 37–86 (subset commonly encountered in bio structures)
    "RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I","XE",
    "CS","BA","LA","CE","PR","ND","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF","TA",
    "W","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN"
}

# Residue names treated as "standard" (kept as ATOM); everything else -> HETATM
STD_RESIDUES = {
    # amino acids (incl. common variants)
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","HID","HIE","HIP","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","SEC","PYL",
    # nucleotides (common forms)
    "A","C","G","T","U","DA","DC","DG","DT","DU","ADE","CYT","GUA","THY","URA",
    "AMP","CMP","GMP","UMP","TMP","ADP","CDP","GDP","UDP","ATP","CTP","GTP","UTP",
}

# Covalent radii (Å) for common elements.
COV_RAD = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39,
    "B": 0.85, "SI": 1.11, "SE": 1.20,
    # Common metals (approximate, single-bond covalent radii)
    "NA": 1.66, "MG": 1.41, "AL": 1.21, "K": 2.03, "CA": 1.76,
    "MN": 1.39, "FE": 1.32, "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22,
    "SR": 1.95, "CD": 1.44, "HG": 1.44, "PT": 1.36, "AU": 1.36, "AG": 1.45,
}

# Bonding thresholds
BOND_SCALE = 1.00
BOND_SLACK = 0.45     # extra Å to tolerate coordinates/rounding
MIN_DIST   = 0.40     # avoid self/near-zero distances

# --- Mapping from PDBQT atom types → PDB element symbols (right-justified 2-chars) ---

def pdbqt_type_to_element(pdbqt_type: str) -> str:
    """
    Convert AutoDock/ADFR PDBQT atom type (e.g., 'A','OA','HD','Cl','Zn') to
    a PDB element symbol (2 chars, right-justified).
    """
    m = {
        # carbons
        "C":"C", "A":"C",                   # A = aromatic carbon nucleus
        # nitrogens
        "N":"N", "NA":"N",
        # oxygens
        "O":"O", "OA":"O", "OW":"O",
        # sulfurs
        "S":"S", "SA":"S",
        # hydrogens
        "H":"H", "HD":"H", "HS":"H",
        # halogens (case variants)
        "F":"F", "CL":"Cl", "Cl":"Cl", "BR":"Br", "Br":"Br", "I":"I",
        # others
        "P":"P", "SI":"Si", "Si":"Si", "SE":"Se", "Se":"Se", "B":"B",
        # metals (allow both cases)
        "MG":"Mg", "Mg":"Mg", "MN":"Mn","Mn":"Mn", "FE":"Fe","Fe":"Fe",
        "CO":"Co","Co":"Co", "NI":"Ni","Ni":"Ni", "CU":"Cu","Cu":"Cu",
        "ZN":"Zn","Zn":"Zn", "CA":"Ca","Ca":"Ca", "NA":"Na","Na":"Na", "K":"K",
        "SR":"Sr","Sr":"Sr", "CD":"Cd","Cd":"Cd", "HG":"Hg","Hg":"Hg",
        "PT":"Pt","Pt":"Pt", "AU":"Au","Au":"Au", "AG":"Ag","Ag":"Ag",
    }
    key = pdbqt_type.strip()
    if key in m:
        el = m[key]
    else:
        # Normalize capitalization and retry
        norm = (key[:1].upper() + key[1:].lower()) if key else ""
        el = m.get(norm, "")
        if not el:
            # fallback to first letter (e.g., "Xx" → "X")
            el = (key[:1].upper() if key else "")
    return el.rjust(2) if el else "  "

# --- Utilities ---

def is_standard_resname(resname: str) -> bool:
    return resname.strip().upper() in STD_RESIDUES

def parse_float(s: str, default: float = 0.0) -> float:
    try:
        return float(s.strip())
    except Exception:
        return default

def fmt_pdb_atom(record: str, serial: int, atom_name: str, resname: str, chain: str,
                 resseq: int, x: float, y: float, z: float, occ: float, temp: float,
                 element: str) -> str:
    atom_name = atom_name[:4].ljust(4)
    resname = resname[:3].rjust(3)
    chain = (chain or "A")[:1]
    return (f"{record:<6s}{serial:>5d}  "
            f"{atom_name}{resname} {chain}"
            f"{resseq:>4d}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{occ:>6.2f}{temp:>6.2f}          "
            f"{element:>2s}\n")

def fmt_conect(i_serial: int, bonded_serials: list[int]) -> list[str]:
    """
    Build one or more CONECT lines for i_serial, grouping at most 4 neighbors per line.
    """
    lines = []
    chunk = []
    for j in bonded_serials:
        chunk.append(j)
        if len(chunk) == 4:
            lines.append("CONECT" + f"{i_serial:>5d}" + "".join(f"{k:>5d}" for k in chunk) + "\n")
            chunk = []
    if chunk:
        lines.append("CONECT" + f"{i_serial:>5d}" + "".join(f"{k:>5d}" for k in chunk) + "\n")
    return lines

def dist2(a, b) -> float:
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return dx*dx + dy*dy + dz*dz

# --- ARG hydrogen name fixer ---

def fix_arg_h_name(atom_name: str) -> str:
    """
    Normalize ARG guanidinium H names to PDB standard:
      1HH1→HH11, 2HH1→HH12, 1HH2→HH21, 2HH2→HH22
    Idempotent for already-correct names.
    """
    name = atom_name.strip().upper()
    mapping = {
        "1HH1": "HH11",
        "2HH1": "HH12",
        "1HH2": "HH21",
        "2HH2": "HH22",
        # standard forms
        "HH11": "HH11",
        "HH12": "HH12",
        "HH21": "HH21",
        "HH22": "HH22",
    }
    return mapping.get(name, atom_name)

# --- Main conversion ---

def convert_pdbqt_to_pdb_with_conect(inp_pdbqt: str, out_pdb: str,
                                     bond_scale: float = BOND_SCALE,
                                     bond_slack: float = BOND_SLACK) -> None:
    atoms = []  # list of dicts with parsed atom data
    serial_next = 1

    # 1) Read PDBQT atoms and map to elements
    with open(inp_pdbqt, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            rec = line[:6].strip().upper()
            if rec in ("ATOM", "HETATM"):
                # fixed-column parse + fallbacks
                serial_str = line[6:11].strip()
                try:
                    serial = int(serial_str) if serial_str else serial_next
                except Exception:
                    serial = serial_next
                atom_name = line[12:16].strip() or "C"
                resname = line[17:20].strip() or "LIG"
                chain = line[21].strip() or "A"
                resseq_str = line[22:26].strip()
                try:
                    resseq = int(resseq_str) if resseq_str else 1
                except Exception:
                    resseq = 1

                # --- NEW: rename non-standard ARG H names instead of dropping ---
                if resname.upper() == "ARG":
                    atom_name = fix_arg_h_name(atom_name)

                x = parse_float(line[30:38], 0.0)
                y = parse_float(line[38:46], 0.0)
                z = parse_float(line[46:54], 0.0)
                occ = parse_float(line[54:60], 1.00)
                temp = parse_float(line[60:66], 0.00)

                # PDBQT atom type (last columns or last token)
                pdbqt_field = line[76:80].strip()
                if not pdbqt_field:
                    toks = line.split()
                    pdbqt_field = toks[-1] if toks else ""
                element = pdbqt_type_to_element(pdbqt_field)  # -> 2-char, right-justified

                record = "ATOM" if is_standard_resname(resname) else "HETATM"

                atoms.append({
                    "record": record,
                    "serial": serial,
                    "name": atom_name,
                    "resname": resname,
                    "chain": chain,
                    "resseq": resseq,
                    "xyz": (x, y, z),
                    "occ": occ,
                    "temp": temp,
                    "element": element.strip() or "C",  # store bare element ('C','O',...)
                })
                serial_next = max(serial_next, serial + 1)
            # ignore PDBQT control records (ROOT, BRANCH, etc.)

    if not atoms:
        raise ValueError("No ATOM/HETATM records found in the PDBQT.")

    # 2) Build bonds only where at least one endpoint is HETATM
    elem_upper = [a["element"].upper() for a in atoms]
    coords = [a["xyz"] for a in atoms]
    is_het = [a["record"] == "HETATM" for a in atoms]

    bonds = set()  # pairs of serials (min, max)
    n = len(atoms)

    # Precompute per-atom covalent radius
    radii = []
    for e in elem_upper:
        r = COV_RAD.get(e, COV_RAD.get(e.upper(), 0.77))  # default ~C
        radii.append(r)

    for i in range(n):
        # If neither atom in a pair is HETATM, we skip; so only check when i is HETATM,
        # to halve the work and guarantee the condition.
        if not is_het[i]:
            continue
        xi = coords[i]; ri = radii[i]
        for j in range(n):
            if i == j:
                continue
            xj = coords[j]; rj = radii[j]
            d2 = dist2(xi, xj)
            if d2 < MIN_DIST * MIN_DIST:
                continue
            # Distance cutoff
            cutoff = bond_scale * (ri + rj) + bond_slack
            if d2 <= cutoff * cutoff:
                a = atoms[i]["serial"]; b = atoms[j]["serial"]
                if a != b:
                    bonds.add((a, b) if a < b else (b, a))

    # 3) Write PDB with CONECT for bonds involving HETATM atoms
    with open(out_pdb, "w", encoding="utf-8") as out:
        # a) ATOM/HETATM
        for a in atoms:
            out.write(fmt_pdb_atom(
                a["record"], a["serial"], a["name"], a["resname"], a["chain"],
                a["resseq"], a["xyz"][0], a["xyz"][1], a["xyz"][2],
                a["occ"], a["temp"], a["element"].rjust(2)
            ))

        # b) CONECT (group neighbors per source atom)
        adj = {}
        for u, v in bonds:
            adj.setdefault(u, []).append(v)
            adj.setdefault(v, []).append(u)

        for a in atoms:
            s = a["serial"]
            neigh = adj.get(s, [])
            if not neigh:
                continue
            neigh_sorted = sorted(neigh)
            for line in fmt_conect(s, neigh_sorted):
                out.write(line)

        out.write("END\n")

# --- Merging helper (ligand + receptor PDBQT → merged PDBQT, then to PDB) ---

def parse_atom_line(line: str) -> Dict:
    """Parse an ATOM/HETATM line (PDB/PDBQT-like). Returns a dict with fields."""
    rec = line[0:6].strip()
    if rec not in ("ATOM", "HETATM"):
        raise ValueError("Not an atom line")

    def safe_slice(s: str, start: int, end: int) -> str:
        return s[start:end] if len(s) >= end else (s[start:] if len(s) > start else "")

    serial_s = safe_slice(line, 6, 11).strip()
    serial = int(serial_s) if serial_s.isdigit() else 0
    name  = safe_slice(line, 12, 16)
    alt   = safe_slice(line, 16, 17)
    resn  = safe_slice(line, 17, 20)
    chain = safe_slice(line, 21, 22)
    resi_s = safe_slice(line, 22, 26).strip()
    resi = int(resi_s) if re.match(r"^-?\d+$", resi_s) else 1
    icode = safe_slice(line, 26, 27)
    x_s = safe_slice(line, 30, 38).strip()
    y_s = safe_slice(line, 38, 46).strip()
    z_s = safe_slice(line, 46, 54).strip()
    occ_s = safe_slice(line, 54, 60).strip()
    b_s   = safe_slice(line, 60, 66).strip()
    elem  = safe_slice(line, 76, 78).strip()
    tail  = safe_slice(line, 66, len(line)).rstrip("\n")

    def _f(v, d):
        try:
            return float(v)
        except Exception:
            return d
    x = _f(x_s, 0.0)
    y = _f(y_s, 0.0)
    z = _f(z_s, 0.0)
    occ = _f(occ_s, 1.00)
    bfac= _f(b_s, 0.00)

    return {
        "record": rec,
        "serial": serial,
        "name": name,
        "altloc": alt,
        "resname": resn,
        "chain": chain,
        "resseq": resi,
        "icode": icode,
        "x": x, "y": y, "z": z,
        "occ": occ, "bfac": bfac,
        "element": elem,
        "tail": tail,
    }

def format_atom_name(name: str) -> str:
    nm = name.strip()
    return f"{nm:<4s}"[:4]

def format_atom_line(a: Dict) -> str:
    """Format atom dict back to PDB-like line, preserving tail (cols ≥67)."""
    elem = a.get("element", "").strip()
    if not elem:
        m = re.match(r"([A-Za-z]{1,2})", a["name"].strip())
        elem = (m.group(1) if m else "")
    elem = elem.upper()[:2]

    line = (
        f"{a['record']:<6}{a['serial']:>5d} "
        f"{format_atom_name(a['name'])}"
        f"{a['altloc'][:1]:1}"
        f"{a['resname'][:3]:>3s} "
        f"{a['chain'][:1]:1}"
        f"{a['resseq']:>4d}{a['icode'][:1]:1}   "
        f"{a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
        f"{a['occ']:>6.2f}{a['bfac']:>6.2f}          "
        f"{elem:>2s}"
    )
    if a.get("tail"):
        t = a["tail"].rstrip()
        if t and not t.startswith(" "):
            line += " "
        line += t
    return line + "\n"

def read_pdbqt_atoms(path: str, origin: str) -> List[Dict]:
    atoms: List[Dict] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    a = parse_atom_line(line)
                    a["origin"] = origin
                    # normalize ARG names here too (if the merge path is used)
                    if a["resname"].strip().upper() == "ARG":
                        a["name"] = fix_arg_h_name(a["name"])
                    atoms.append(a)
                except Exception:
                    continue
    return atoms

def write_pdbqt_atoms(path: str, atoms: List[Dict]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for a in atoms:
            out.write(format_atom_line(a))

def normalize_ligand_residue_and_chain(lig: List[Dict]) -> Tuple[str, str]:
    if not lig:
        return ("", "")
    first = lig[0]
    ref_resn = first["resname"].strip() or "LIG"
    ref_chain = first["chain"].strip() or "A"
    for a in lig:
        a["resname"] = f"{ref_resn:>3s}"[:3]
        a["chain"] = f"{ref_chain:1s}"[:1]
    return (ref_resn, ref_chain)

def deduplicate_receptor_atoms(rec: List[Dict], lig: List[Dict]) -> List[Dict]:
    lig_keys = {(a["name"].strip(), a["resname"].strip(), a["chain"].strip()) for a in lig}
    filtered: List[Dict] = []
    for a in rec:
        key = (a["name"].strip(), a["resname"].strip(), a["chain"].strip())
        if key in lig_keys:
            continue
        filtered.append(a)
    return filtered

def looks_like_element_only(name: str) -> str | None:
    n = name.strip()
    m = re.fullmatch(r"([A-Za-z]{1,2})", n)
    if not m:
        return None
    cand = m.group(1).upper()
    return cand if cand in ELEMENTS else None

def rename_weak_ligand_atom_names(lig: List[Dict], existing: List[Dict]) -> None:
    """For ligand atoms whose name is just the element, rename to Element+index within each residue."""
    used: Dict[Tuple[str,int,str], set] = {}
    def reskey(a: Dict) -> Tuple[str,int,str]:
        return (a["resname"].strip(), int(a["resseq"]), a["chain"].strip())

    for a in existing:
        rk = reskey(a)
        used.setdefault(rk, set()).add(a["name"].strip().upper())

    for a in lig:
        rk = reskey(a)
        used.setdefault(rk, set())
        elem = looks_like_element_only(a["name"])
        if elem is None:
            nm = a["name"].strip().upper()
            if nm in used[rk]:
                base = (a["element"].strip().upper() or nm[:2] or "X")
                idx = 1
                while True:
                    cand = f"{base}{idx}"[:4]
                    if cand not in used[rk]:
                        a["name"] = cand
                        a["element"] = base[:2]
                        used[rk].add(cand)
                        break
                    idx += 1
            else:
                used[rk].add(nm)
            continue

        base = elem
        idx = 1
        while True:
            cand = f"{base}{idx}"[:4]
            if cand not in used[rk]:
                a["name"] = cand
                a["element"] = base[:2]
                used[rk].add(cand)
                break
            idx += 1

def _which_openbabel() -> str | None:
    for exe in ("obabel", "babel"):
        p = shutil.which(exe)
        if p:
            return p
    return None

def convert_with_openbabel(in_pdbqt: str, out_prefix: str) -> None:
    exe = _which_openbabel()
    if not exe:
        print("[WARN] Open Babel (obabel/babel) not found on PATH. Skipping conversion.")
        return

    mol2_path = f"{out_prefix}.mol2"
    pdb_path  = f"{out_prefix}.pdb"

    # PDBQT → MOL2
    cmd1 = [exe, "-ipdbqt", in_pdbqt, "-omol2", "-O", mol2_path]
    print("[RUN] ", " ".join(cmd1))
    subprocess.run(cmd1, check=True)

    # MOL2 → PDB
    cmd2 = [exe, "-imol2", mol2_path, "-opdb", "-O", pdb_path]
    print("[RUN] ", " ".join(cmd2))
    subprocess.run(cmd2, check=True)

def run(lig_path: str, rec_path: str, out_path: str) -> None:
    if not os.path.isfile(lig_path):
        raise FileNotFoundError(f"Ligand not found: {lig_path}")
    if not os.path.isfile(rec_path):
        raise FileNotFoundError(f"Receptor not found: {rec_path}")

    lig = read_pdbqt_atoms(lig_path, origin="ligand")
    rec = read_pdbqt_atoms(rec_path, origin="receptor")

    if not lig:
        raise RuntimeError("Ligand file contains no ATOM/HETATM records.")

    # Step 1: normalize ligand resName/chainID
    ref_resn, ref_chain = normalize_ligand_residue_and_chain(lig)
    lig_norm_path = lig_path+'_norm'
    write_pdbqt_atoms(lig_norm_path, lig)
    print(f"[OK] Wrote normalized ligand → {lig_norm_path} (resName={ref_resn}, chain={ref_chain})")

    # Step 2: de-duplicate receptor against ligand
    rec_filtered = deduplicate_receptor_atoms(rec, lig)
    print(f"[OK] Receptor atoms: {len(rec)} → {len(rec_filtered)} after duplicate removal")

    # Step 3: rename weak ligand atom names to avoid clashes
    rename_weak_ligand_atom_names(lig, existing=rec_filtered)

    # Merge (receptor first, then ligand)
    merged = list(rec_filtered) + list(lig)

    # Re-assign serials 1..N for cleanliness
    for i, a in enumerate(merged, start=1):
        a["serial"] = i

    merged_pdbqt = out_path
    write_pdbqt_atoms(merged_pdbqt, merged)
    print(f"[OK] Wrote merged PDBQT → {merged_pdbqt}  (total atoms: {len(merged)})")

def combine_pdbqt_convert_pdb(lig, rec, out, out_pdb):
    run(lig, rec, out)
    convert_pdbqt_to_pdb_with_conect(out, out_pdb)

