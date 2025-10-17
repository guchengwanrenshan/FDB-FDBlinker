#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ftlk.py – Flexible two-ligand covalent docking launcher (Py3 main)

This script:
  • Uses Open Babel 3 (pybel) from the active conda env for all OB work in Py3.
  • Calls legacy ADFR/MGLTools utilities under ADFR’s own Python-2 runtimes:
      - pythonsh (Py2) for prepare_* and MolKit cleanup simulation.
      - python2.7 (Py2) for runAGFR.py and runADFR.py.
  • Forces ADFR’s Open Babel 2.4.1 for the Py2 subprocesses only.

Environment overrides (optional):
  ADFR_HOME                default: /home/zl/app/ADFRsuite_x86_64Linux_1.0
  ADFR_PY                  default: $ADFR_HOME/bin/python2.7
  ADFR_PY_MGL              default: $ADFR_HOME/bin/pythonsh
  CCSB_ROOT                default: $ADFR_HOME/CCSBpckgs
  PYTHON_EXE               default: sys.executable (current Py3)
  FTLK_HELPER_PY           default: $ADFR_HOME/CCSBpckgs/ADFR/bin/partial_minimize_ob3.py
  FTLK_HELPER_PYTHON       default: PYTHON_EXE (Py3)
  FTLK_HELPER_DISABLE=1    disable the optional OB3 minimizer helper

Usage (example):
  conda activate ob3
  export ADFR_HOME="/home/zl/app/ADFRsuite_x86_64Linux_1.0"
  export ADFR_PY="$ADFR_HOME/bin/python2.7"
  export ADFR_PY_MGL="$ADFR_HOME/bin/pythonsh"
  export CCSB_ROOT="$ADFR_HOME/CCSBpckgs"
  ./ftlk --complex "/home/zl/app/uploads/0.pdb" \
         --ligand  "/home/zl/app/uploads/Toluene.pdbqt" \
         --atomsA 4477 4474 4478 --atomsB 12 \
         --jobName colvant_t.pdb --outputPath /home/zl/app/uploads \
         --maxEvals 400000 --nbRuns 8
"""

from __future__ import annotations
import os, sys, argparse, subprocess, tempfile, shutil, traceback, re, glob, math, itertools, textwrap, tarfile
from pathlib import Path
from pdbqt_tools import combine_pdbqt_convert_pdb

# ======== Interpreter + path configuration ========

ADFR_HOME   = os.environ.get("ADFR_HOME", "/home/zl/app/ADFRsuite_x86_64Linux_1.0")
ADFR_PY     = os.environ.get("ADFR_PY",     f"{ADFR_HOME}/bin/python2.7")   # AGFR/ADFR
ADFR_PY_MGL = os.environ.get("ADFR_PY_MGL", f"{ADFR_HOME}/bin/pythonsh")    # prepare_* + MolKit helper
CCSB_ROOT   = os.environ.get("CCSB_ROOT",   f"{ADFR_HOME}/CCSBpckgs")       # MolKit/ADT tree

PYTHON_EXE    = os.environ.get("PYTHON_EXE", sys.executable)                # current Py3
HELPER_PYTHON = os.environ.get("FTLK_HELPER_PYTHON", PYTHON_EXE)
HELPER_PY     = os.environ.get("FTLK_HELPER_PY", f"{ADFR_HOME}/CCSBpckgs/ADFR/bin/partial_minimize_ob3.py")
HELPER_DISABLE= os.environ.get("FTLK_HELPER_DISABLE", "")

print("[ftlk] Py3 main:", PYTHON_EXE)
print("[ftlk] Py2 ADFR :", ADFR_PY)
print("[ftlk] Py2 MGL  :", ADFR_PY_MGL)

PKG_BIN = os.path.join(ADFR_HOME, "CCSBpckgs", "ADFR", "bin")
_RX_SKIP = re.compile(r"^(ROOT|ENDROOT|BRANCH|ENDBRANCH|TORSDOF|REMARK)")

# ======== Prefer Open Babel 3 from the active conda env (for Py3) ========

def _prefer_conda_openbabel():
    cp = os.environ.get("CONDA_PREFIX", "")
    if not cp: return
    cp = Path(cp)
    # Data
    for d in (cp/"share/openbabel", cp/"Library/share/openbabel"):
        if d.is_dir(): os.environ["BABEL_DATADIR"] = str(d); break
    # Plugins
    for d in (cp/"lib/openbabel/3.1", cp/"lib/openbabel/3.1.1",
              cp/"Library/lib/openbabel/3.1", cp/"Library/lib/openbabel/3.1.1"):
        if d.is_dir(): os.environ["OB_PLUGIN_PATH"] = str(d); break
    # PATH precedence
    parts = os.environ.get("PATH","").split(os.pathsep)
    for b in map(str, (cp/"bin", cp/"Library/bin", cp/"Scripts")):
        if Path(b).is_dir() and b not in parts: parts.insert(0, b)
    os.environ["PATH"] = os.pathsep.join(parts)

_prefer_conda_openbabel()

try:
    from openbabel import openbabel as ob
    from openbabel import pybel
except Exception:
    try:
        import pybel; from openbabel import openbabel as ob  # type: ignore
    except Exception as e:
        print("[ftlk] FATAL: Open Babel (pybel) not found in this Python 3 environment.")
        raise

def _get_atomic_num(symbol: str) -> int:
    symbol = symbol.strip().title()
    try:
        if hasattr(ob, "etab"): return ob.etab.GetAtomicNum(symbol)
        if hasattr(ob, "OBElementTable"): return ob.OBElementTable().GetAtomicNum(symbol)
    except Exception:
        pass
    PT = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,"Na":11,"Mg":12,"Al":13,"Si":14,
          "P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,
          "Ni":28,"Cu":29,"Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,"Rb":37,"Sr":38,"Y":39,"Zr":40,
          "Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,"In":49,"Sn":50,"Sb":51,"Te":52,"I":53,"Xe":54}
    return PT.get(symbol, 6)

# ======== ADFR’s OB-2.4 environment for Py2 subprocesses (AGFR/ADFR) ========
'''
def _adfr_ob_env(base_env=None):
    """
    Env for ADFR's Python2 subprocesses (runAGFR/runADFR) using Open Babel 2.4.1.
    Key: OB-2.4 uses BABEL_LIBDIR (NOT OB_PLUGIN_PATH).
    """
    from pathlib import Path
    env = dict(base_env or os.environ)
    A = ADFR_HOME
    root = Path(A)

    # Confirmed on your install:
    babel_datadir = root / "share" / "openbabel" / "2.4.1"
    babel_libdir  = root / "lib"   / "openbabel" / "2.4.1"
    adfr_bin      = root / "bin"
    adfr_lib      = root / "lib"

    # (Some bundles also ship 3rd-party libs separately)
    thirdparty_lib = root / "ThirdPartyPacks" / "lib"

    # Required for OB-2.4:
    env["BABEL_DATADIR"] = str(babel_datadir)
    env["BABEL_LIBDIR"]  = str(babel_libdir)   # <-- critical for 2.4

    # Make ADFR bins/libs (and 3rd-party libs) take precedence
    env["PATH"] = str(adfr_bin) + os.pathsep + env.get("PATH", "")
    ld_parts = [str(adfr_lib)]
    if thirdparty_lib.is_dir():
        ld_parts.append(str(thirdparty_lib))
    if env.get("LD_LIBRARY_PATH"):
        ld_parts.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = os.pathsep.join(ld_parts)

    # ADFR Python packages
    env["PYTHONPATH"] = f"{A}/CCSBpckgs" + os.pathsep + env.get("PYTHONPATH", "")

    # Scrub OB-3 variables that can confuse OB-2.4
    for k in list(env.keys()):
        ku = k.upper()
        if ku in ("OB_PLUGIN_PATH",):  # OB3 var – drop it for Py2 tools
            del env[k]
        # Keep our OB-2.4 vars; they start with BABEL_
    return env
'''
def _adfr_ob_env(base_env=None):
    env = dict(base_env or os.environ)
    A = ADFR_HOME
    env["BABEL_DATADIR"]   = f"{A}/share/openbabel/2.4.1"
    env["BABEL_LIBDIR"]    = f"{A}/lib/openbabel/2.4.1"
    env["PATH"]            = f"{A}/bin" + os.pathsep + env.get("PATH","")
    env["LD_LIBRARY_PATH"] = f"{A}/lib" + os.pathsep + env.get("LD_LIBRARY_PATH","")
    env["PYTHONPATH"]      = f"{A}/CCSBpckgs" + os.pathsep + env.get("PYTHONPATH","")

    # 👉 Tell ADFR exactly where pythonsh is:
    env["PYTHONSH"] = f"{A}/bin/pythonsh"
    env["MGLPYTHON"] = env["PYTHONSH"]  # harmless; some scripts check this

    # scrub OB3 vars
    for k in list(env):
        if k.upper().startswith("OB_"):
            del env[k]
    return env

# ======== MolKit/ADT helper (Py2 via pythonsh) to get post-cleanup root index ========

def _mgl_find_root_index(pdb_path: str, anchor_name: str, anchor_resname: str) -> int:
    code = textwrap.dedent("""
        import sys, json
        from MolKit import Read
        from AutoDockTools.MoleculePreparation import AutoDockMoleculePreparation
        pdb, aname, rname = sys.argv[1], sys.argv[2], sys.argv[3]
        mol = Read(pdb)[0]
        AutoDockMoleculePreparation(mol, mode='automatic', repairs='',
                                    charges_to_add='gasteiger', cleanup='nphs_lps',
                                    outputfilename=None, debug=False, version=4)
        idx = None
        for i, at in enumerate(mol.allAtoms):
            res = getattr(getattr(at, 'parent', None), 'type', '')
            if at.name == aname and res == rname:
                idx = i; break
        print(json.dumps({'ok': idx is not None, 'idx': idx}))
    """)
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
        hp = f.name; f.write(code)
    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = CCSB_ROOT + os.pathsep + env.get("PYTHONPATH","")
        out = subprocess.check_output([ADFR_PY_MGL, hp, pdb_path, anchor_name, anchor_resname], env=env, text=True)
        import json
        data = json.loads(out.strip().splitlines()[-1])
        if not data.get("ok"): raise RuntimeError("MolKit helper could not locate anchor")
        return int(data["idx"])
    finally:
        try: os.remove(hp)
        except OSError: pass

# ======== File/build helpers ========

def create_complex_receptor_file(original_pdb: str, final_receptor_path: str) -> None:
    print('[ftlk] Preparing complex file...')
    workdir = os.path.dirname(final_receptor_path)
    protein_pdb_path = os.path.join(workdir, 'protein_original.pdb')
    ligand_pdb_path  = os.path.join(workdir, 'ligandA_original.pdb')

    # Split ATOM (protein) and HETATM (LigandA) from the complex PDB
    with open(original_pdb, 'r') as fin, \
         open(protein_pdb_path, 'w') as fp, \
         open(ligand_pdb_path,  'w') as fl:
        for line in fin:
            if line.startswith('ATOM'): fp.write(line)
            elif line.startswith('HETATM'): fl.write(line)
            elif line.startswith('TER'): fp.write(line)

    if not os.path.getsize(ligand_pdb_path) > 0:
        raise RuntimeError('HETATM records for LigandA not found in %s.' % original_pdb)

    mgl_util = os.path.join(CCSB_ROOT, 'AutoDockTools', 'Utilities24')
    prep_rec_script = os.path.join(mgl_util, 'prepare_receptor4.py')
    prep_lig_script = os.path.join(mgl_util, 'prepare_ligand4.py')
    protein_qt_path = protein_pdb_path.replace('.pdb', '.pdbqt')
    ligand_qt_path  = ligand_pdb_path.replace('.pdb',  '.pdbqt')

    env = os.environ.copy()
    env["PYTHONPATH"] = CCSB_ROOT + os.pathsep + env.get("PYTHONPATH","")

    # Prepare receptor (Py2 pythonsh)
    subprocess.check_call([ADFR_PY_MGL, prep_rec_script,
                           '-r', os.path.basename(protein_pdb_path),
                           '-A', 'checkhydrogens', '-U', 'nphs_lps_waters',
                           '-o', os.path.basename(protein_qt_path)], cwd=workdir, env=env)

    # Prepare LigandA (Py2 pythonsh)
    subprocess.check_call([ADFR_PY_MGL, prep_lig_script,
                           '-l', os.path.basename(ligand_pdb_path),
                           '-A', 'checkhydrogens',
                           '-o', os.path.basename(ligand_qt_path)], cwd=workdir, env=env)

    # Merge into one receptor PDBQT (renumbering)
    final_lines, total = [], 0
    last_prot = None
    with open(protein_qt_path, 'r') as f:
        for line in f:
            if _RX_SKIP.match(line) or line.startswith('TER'): continue
            if line.startswith('ATOM'):
                total += 1
                final_lines.append('%s%5d%s' % (line[:6], total, line[11:]))
                last_prot = line
    if last_prot:
        res_name = last_prot[17:20]; chain_id = last_prot[21:22]; res_seq = last_prot[22:26]
        final_lines.append('TER   %5d      %s %s%s\n' % (total + 1, res_name, chain_id, res_seq))
        total += 1
    with open(ligand_qt_path, 'r') as f:
        for line in f:
            if _RX_SKIP.match(line): continue
            if line.startswith('HETATM'):
                total += 1
                final_lines.append('ATOM  %5d%s' % (total, line[11:]))

    with open(final_receptor_path, 'w') as fout:
        fout.write(''.join(final_lines))
    print('[ftlk] Final receptor file:', final_receptor_path, f'({total} atoms)')

def find_new_anchor_serials(original_pdb_path: str, new_pdbqt_path: str, original_serial_list: list[int]) -> list[int]:
    new_serials = []
    with open(new_pdbqt_path, 'r') as f_new: lines = f_new.readlines()
    for orig_serial in original_serial_list:
        atom_name = resname = resseq = None
        with open(original_pdb_path, 'r') as f_orig:
            for line in f_orig:
                if line.startswith('HETATM'):
                    try:
                        if int(line[6:11]) == orig_serial:
                            atom_name = line[12:16].strip()
                            resname   = line[17:20].strip()
                            resseq    = line[22:26].strip()
                            break
                    except Exception:
                        continue
        if atom_name is None:
            raise RuntimeError("HETATM serial %d not found in %s." % (orig_serial, original_pdb_path))
        found = None
        for line in lines:
            if line.startswith(('ATOM','HETATM')):
                try:
                    if (line[12:16].strip()==atom_name and line[17:20].strip()==resname and line[22:26].strip()==resseq):
                        found = int(line[6:11]); break
                except Exception:
                    continue
        if found is None:
            raise RuntimeError("Atom '%s' %s%s not found in %s." % (atom_name, resname, resseq, new_pdbqt_path))
        new_serials.append(found)
    return new_serials

def get_atom_coords_by_serial(pdbqt_path: str, serial_number: int) -> tuple[float,float,float]:
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM','HETATM')) and int(line[6:11]) == serial_number:
                return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    raise ValueError('Atom %d not found in %s.' % (serial_number, pdbqt_path))

def find_residue_info(pdbqt_path: str, atom_serial: int) -> tuple[str,str,str]:
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM','HETATM')) and int(line[6:11]) == atom_serial:
                chain_id = line[21:22].strip() or ' '
                res_name = line[17:20].strip()
                res_seq  = line[22:26].strip()
                return chain_id, res_name, res_seq
    raise ValueError('Atom %d not found in %s.' % (atom_serial, pdbqt_path))

def _get_atom_info_manually(pdbqt_path: str, serial_to_find: int) -> dict:
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM','HETATM')):
                try:
                    if int(line[6:11]) == serial_to_find:
                        return {
                            'name': line[12:16].strip(),
                            'res_name': line[17:20].strip(),
                            'chain': line[21:22].strip(),
                            'res_num': int(line[22:26].strip())
                        }
                except Exception:
                    continue
    raise ValueError("Atom serial %d not found in %s." % (serial_to_find, pdbqt_path))

def _fix_sidechain_info_in_pdb(file_path: str, sidechain_atom_names: list[str], sidechain_atom_infos: list[dict]) -> None:
    with open(file_path, 'r') as f: lines = f.readlines()
    num_side = len(sidechain_atom_names)
    atom_idx = [i for i, ln in enumerate(lines) if ln.startswith(('ATOM','HETATM'))]
    targets  = [i for i in atom_idx if not lines[i][12:16].strip()]
    if len(targets) != num_side:
        raise RuntimeError("Sidechain correction mismatch: expected %d blanks, found %d." % (num_side, len(targets)))
    for i in range(num_side):
        li = targets[i]; line = lines[li]
        aname = sidechain_atom_names[i]
        line = line[:12] + (" " + aname).ljust(4) + line[16:]
        info = sidechain_atom_infos[i]
        res_segment = "{:>3s} {:1s}{:>4s} ".format(info['res_name'][:3].strip(), info['chain'][:1].strip(), str(info['res_num']))
        lines[li] = line[:17] + res_segment + line[27:]
    for i, ln in enumerate(lines):
        if ln.startswith('ATOM  '): lines[i] = "HETATM" + ln[6:]
    with open(file_path, 'w') as f: f.writelines(lines)

def _restore_atom_record_type_in_pdbqt(pdbqt_path: str, num_atoms_to_change: int = 3) -> None:
    with open(pdbqt_path, 'r') as f: lines = f.readlines()
    changed = 0
    for i, ln in enumerate(lines):
        if changed >= num_atoms_to_change: break
        if ln.startswith('HETATM'):
            lines[i] = "ATOM  " + ln[6:]; changed += 1
    with open(pdbqt_path, 'w') as f: f.writelines(lines)

def _swap_root_atoms_and_fix_branch(pdbqt_path: str) -> None:
    try:
        with open(pdbqt_path, 'r') as f: lines = f.readlines()
        root = next((i for i, ln in enumerate(lines) if ln.startswith('ROOT')), -1)
        if root == -1: return
        a1, a2 = root+1, root+2
        if not (lines[a1].startswith(('ATOM','HETATM')) and lines[a2].startswith(('ATOM','HETATM'))): return
        c1, c2 = lines[a1][11:], lines[a2][11:]
        lines[a1] = lines[a1][:11] + c2; lines[a2] = lines[a2][:11] + c1
        endroot = next((i for i in range(root, len(lines)) if lines[i].startswith('ENDROOT')), -1)
        if endroot != -1:
            for i in range(endroot+1, len(lines)):
                if lines[i].startswith('BRANCH'):
                    try:
                        if int(lines[i][6:10].strip()) == 1:
                            lines[i] = "BRANCH%4d%s" % (2, lines[i][10:]); break
                    except Exception:
                        continue
        with open(pdbqt_path, 'w') as f: f.writelines(lines)
    except Exception as e:
        print("[ftlk] WARNING: swap/fix failed:", e)

# ======== Optional OB3 minimization helper (Py3) ========

def _run_ob3_helper(mol2_in: str, pdb_out: str, movable_indices: list[int], num_new: int) -> bool:
    if HELPER_DISABLE:
        print("[ftlk] OB3 helper disabled."); return False
    helper_py = os.path.expanduser(os.path.expandvars(HELPER_PY))
    py3 = os.path.expanduser(os.path.expandvars(HELPER_PYTHON))
    if not Path(helper_py).is_file():
        print(f"[ftlk] Helper script not found: {helper_py}. Skipping."); return False

    movable_csv = ",".join(str(i) for i in movable_indices)
    env = os.environ.copy()

    # Favor the current conda OB3 for the helper
    pyroot = Path(py3).parent
    path = [str(pyroot), str(pyroot/"Library/bin"), str(pyroot/"Scripts")] + \
           [p for p in env.get("PATH","").split(os.pathsep) if "OpenBabel-2.4.1" not in p]
    env["PATH"] = os.pathsep.join([p for p in path if Path(p).is_dir()])
    for cand in (pyroot/"share/openbabel", pyroot/"Library/share/openbabel"):
        if cand.is_dir(): env['BABEL_DATADIR'] = str(cand); break
    for cand in (pyroot/"lib/openbabel/3.1", pyroot/"lib/openbabel/3.1.1",
                 pyroot/"Library/lib/openbabel/3.1", pyroot/"Library/lib/openbabel/3.1.1"):
        if cand.is_dir(): env['OB_PLUGIN_PATH'] = str(cand); break

    res = subprocess.run(
        [py3, helper_py, "--in", mol2_in, "--out", pdb_out,
         "--movable", movable_csv, "--steps", "250", "--ff", "mmff94", "--num-new", str(num_new)],
        env=env, capture_output=True, text=True
    )
    if res.returncode != 0:
        print("[ftlk] OB3 helper stderr:\n", res.stderr or "(empty)")
        print("[ftlk] OB3 helper stdout:\n", res.stdout or "(empty)")
        print(f"[ftlk] WARNING: helper failed with code {res.returncode}.")
        return False
    return True

# ======== Build covalent ligand (OB3 in Py3, prepare_ligand4.py in Py2) ========


def create_merged_ligand_pdbqt(ligand_b_path: str, complex_path: str, atoms_a_serials: list[int],
                               atom_b_serial: int, final_output_path: str, workdir: str, python_exe: str) -> list[int]:
    print('[ftlk] Creating covalent ligand...')
    try:
        ligand_b_mol = next(pybel.readfile("pdbqt", ligand_b_path))
    except StopIteration:
        raise IOError("Open Babel could not read LigandB PDBQT.")

    side = pybel.ob.OBMol()
    map_serial_to_idx, side_names, side_infos = {}, [], []
    for serial in atoms_a_serials:
        info = _get_atom_info_manually(complex_path, serial); side_infos.append(info)
        aname = info['name']; side_names.append(aname)
        def _element_from_pdb_atom_name(aname: str) -> str:
            """Infer element from a PDB-style atom name (e.g., CA, CB, SD, 1HH1).
            Rules:
              - If it begins with Cl/BR/Si (case-insensitive), return Cl/Br/Si.
              - Else return the FIRST alphabetic character uppercased (C/N/O/S/H/P/F/I/B/...).
              - Fallback to 'C'.
            """
            n = (aname or "").strip()
            if len(n) >= 2 and n[:2].upper() in ("CL", "BR"):
                return n[:1].upper() + n[1].lower()   # Cl / Br
            if len(n) >= 2 and n[:2].upper() == "SI":
                return "Si"
            for ch in n:
                if ch.isalpha():
                    return ch.upper()
            return "C"

        # ... inside create_merged_ligand_pdbqt(), replace the old block:
        #   symbol = aname.strip().title()[:2]
        #   if symbol not in ("Cl","Br","Si","Na","Mg","Al","Ca","Fe","Zn","Sn"): symbol = aname[0].upper()
        #   Z = _get_atomic_num(symbol)
        # with:
        symbol = _element_from_pdb_atom_name(aname)
        Z = _get_atomic_num(symbol)
     
        a = side.NewAtom(); a.SetAtomicNum(Z)
        res = side.NewResidue(); res.SetName(info['res_name']); res.SetNum(info['res_num']); res.SetChain(info['chain'])
        res.AddAtom(a); res.SetAtomID(a, str(aname))
        map_serial_to_idx[serial] = a.GetIdx()

    if len(atoms_a_serials) > 1:
        side.AddBond(map_serial_to_idx[atoms_a_serials[1]], map_serial_to_idx[atoms_a_serials[0]], 1)
    if len(atoms_a_serials) > 2:
        side.AddBond(map_serial_to_idx[atoms_a_serials[0]], map_serial_to_idx[atoms_a_serials[2]], 1)

    nB = ligand_b_mol.OBMol.NumAtoms()
    ligand_b_mol.OBMol += side

    attach_b_idx = ligand_b_mol.atoms[atom_b_serial - 1].idx
    attach_a_idx = map_serial_to_idx[atoms_a_serials[1]] + nB
    ligand_b_mol.OBMol.AddBond(attach_b_idx, attach_a_idx, 1)

    for atom in pybel.ob.OBMolAtomIter(ligand_b_mol.OBMol):
        r = atom.GetResidue()
        if r: r.SetName("LIG"); r.SetNum(1); r.SetChain("A")

    mol2_path = os.path.join(workdir, 'merged_ligand_tmp.mol2')
    ligand_b_mol.write("mol2", mol2_path, overwrite=True)

    pdb_path = os.path.join(workdir, 'merged_ligand.pdb')
    movable = [nB + i for i in range(1, len(atoms_a_serials)+1)]
    if not _run_ob3_helper(mol2_path, pdb_path, movable, len(atoms_a_serials)):
        mol_tmp = next(pybel.readfile("mol2", mol2_path))
        mol_tmp.write("pdb", pdb_path, overwrite=True)

    _fix_sidechain_info_in_pdb(pdb_path, side_names, side_infos)
    print('[ftlk] Finding post-cleanup root index (MolKit via pythonsh)...')
    root_idx = _mgl_find_root_index(pdb_path, side_names[2], side_infos[2]['res_name'])

    print('[ftlk] Preparing ligand PDBQT (pythonsh)...')
    mgl_util = os.path.join(CCSB_ROOT, 'AutoDockTools', 'Utilities24')
    prep_lig = os.path.join(mgl_util, 'prepare_ligand4.py')
    env = os.environ.copy(); env["PYTHONPATH"] = CCSB_ROOT + os.pathsep + env.get("PYTHONPATH","")
    subprocess.check_call([ADFR_PY_MGL, prep_lig,
                           '-l', os.path.basename(pdb_path),
                           '-o', os.path.basename(final_output_path),
                           '-A', 'checkhydrogens',
                           '-R', str(root_idx),
                           '-v'], cwd=workdir, env=env)

    _restore_atom_record_type_in_pdbqt(final_output_path, num_atoms_to_change=len(atoms_a_serials))
    _swap_root_atoms_and_fix_branch(final_output_path)
    return [1, 2, 3]

# ======== Postprocessing (clean PDBQT & build complex PDB) ========

def process_model(model_lines: list[str]) -> list[str]:
    processed, atom_map, keep = [], {}, []
    for line in model_lines:
        if line.startswith(('ATOM','HETATM')):
            try:
                if int(line[6:11]) > 3: keep.append(line)
            except Exception:
                continue
    if not keep:
        for line in model_lines:
            if not line.startswith(('ATOM','HETATM','ROOT','ENDROOT','BRANCH','ENDBRANCH','TORSDOF')):
                processed.append(line)
        return processed
    for i, line in enumerate(keep): atom_map[int(line[6:11])] = i+1
    root_placed = False
    for line in model_lines:
        if line.startswith(('MODEL','USER','REMARK')):
            processed.append(line); continue
        if line.startswith(('ATOM','HETATM')):
            try:
                s = int(line[6:11])
                if s in atom_map:
                    nl = "{:6s}{:5d}{}".format(line[:6], atom_map[s], line[11:])
                    if not root_placed:
                        processed += ["ROOT\n", nl, "ENDROOT\n"]; root_placed=True
                    else:
                        processed.append(nl)
            except Exception:
                continue
        elif line.startswith(('BRANCH','ENDBRANCH')):
            try:
                parts = line.split(); kw, s1, s2 = parts[0], int(parts[1]), int(parts[2])
                if s1 in atom_map and s2 in atom_map:
                    processed.append("{:6s}{:4d}{:4d}\n".format(kw, atom_map[s1], atom_map[s2]))
            except Exception:
                continue
        elif line.startswith('ENDMDL'):
            processed.append(line)
    return processed

def run_postprocessing(input_file: str, output_file: str) -> None:
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        chunk = []
        for line in infile:
            if line.startswith('MODEL'):
                if chunk: outfile.writelines(process_model(chunk))
                chunk = [line]
            elif line.startswith('ENDMDL'):
                if chunk:
                    chunk.append(line); outfile.writelines(process_model(chunk)); chunk=[]
            elif chunk:
                chunk.append(line)
            else:
                outfile.write(line)
        if chunk: outfile.writelines(process_model(chunk))

def calculate_ligand_diameter(ligand_pdbqt_path: str) -> float:
    mol = next(pybel.readfile("pdbqt", ligand_pdbqt_path))
    coords = [a.coords for a in mol.atoms]
    if len(coords) < 2: return 0.0
    md2 = 0.0
    for c1, c2 in itertools.combinations(coords, 2):
        d2 = (c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2
        if d2 > md2: md2 = d2
    return math.sqrt(md2)

def extract_tar(tar_path, extract_path):
    """
    Extract a .tar, .tar.gz, .tar.bz2, or .tar.xz archive.
    
    Args:
        tar_path (str): Path to the tar file.
        extract_path (str): Directory to extract to (default: current dir).
    """
    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(path=extract_path)
        print(f"Extracted {tar_path} to {extract_path}")

# ======== Main ========

def main():
    p = argparse.ArgumentParser(description='Two-ligand covalent docking wrapper (AGFR two-pass).',
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--complex', required=True, help='Protein+LigandA complex (PDB)')
    p.add_argument('--ligand',  required=True, help='LigandB to dock (PDBQT)')
    p.add_argument('--atomsA',  nargs=3, type=int, required=True,
                  metavar=('BOND_ATOM_1','BOND_ATOM_2','TORSION_ANCHOR'))
    p.add_argument('--atomsB',  type=int, required=True, metavar=('ATTACH_ATOM'))
    p.add_argument('--box_size', type=float, default=None)
    p.add_argument('--nbRuns',   type=int, default=8)
    p.add_argument('--maxEvals', type=int, default=250000)
    p.add_argument('--seed',     type=int, default=1)
    p.add_argument('--keepTmp', action='store_true')
    p.add_argument('--jobName', type=str, default=None)
    p.add_argument('--outputPath', type=str, default=None)
    args, extra = p.parse_known_args()

    if not args.complex.lower().endswith('.pdb'):
        print("[ftlk] ERROR: --complex must be a PDB file."); sys.exit(1)

    out_dir = args.outputPath
    if out_dir: os.makedirs(out_dir, exist_ok=True); print('[ftlk] Output dir:', os.path.abspath(out_dir))

    tmpdir = tempfile.mkdtemp(prefix='ftlk_'); print('[ftlk] Temp dir:', tmpdir)

    try:
        # 1) Prepare receptor (Py2 pythonsh), merge to PDBQT
        rec_pdbqt = os.path.join(tmpdir, 'complex_asReceptor.pdbqt')
        create_complex_receptor_file(args.complex, rec_pdbqt)

        # 2) Map original LigandA atoms -> new serials in merged receptor
        newA = find_new_anchor_serials(args.complex, rec_pdbqt, args.atomsA)
        print('[ftlk] Translated LigandA serials', args.atomsA, '->', newA)

        # 3) Build covalent ligand (OB3 in Py3), prepare ligand PDBQT with correct root (Py2)
        lig_pdbqt = os.path.join(tmpdir, 'ligandB_with_sidechain.pdbqt')
        newC = create_merged_ligand_pdbqt(args.ligand, rec_pdbqt, newA, args.atomsB, lig_pdbqt, tmpdir, PYTHON_EXE)
        print('[ftlk] ADFR -C serials:', newC)

        # 4) Box
        center = get_atom_coords_by_serial(rec_pdbqt, newA[0])
        if args.box_size is None:
            dia = calculate_ligand_diameter(args.ligand)
            box_dim = dia * 2 + 2.0
        else:
            box_dim = float(args.box_size)
        box = (box_dim, box_dim, box_dim)
        print('[ftlk] Box center: %.3f %.3f %.3f' % center)
        print('[ftlk] Box size  : %.3f %.3f %.3f' % box)

        # 5) AGFR target (Py2) with ADFR’s OB2.4.1 env
        env_py2 = _adfr_ob_env()
        # Optional quick probe so OB2.4 failures are readable
        probe = subprocess.run(
            [ADFR_PY, "-c",
             "import os,openbabel,pybel;"
             "print('OB2.4 probe');"
             "print('BABEL_DATADIR=',os.getenv('BABEL_DATADIR'));"
             "print('BABEL_LIBDIR=',os.getenv('BABEL_LIBDIR'));"
             "c=openbabel.OBConversion();"
             "s=c.GetSupportedInputFormat();"
             "print('formats sample=', s[:80])"],
            env=env_py2, capture_output=True, text=True
        )
        if probe.returncode != 0:
            print("[ftlk] OB2.4 probe stderr:\n", probe.stderr or "(empty)")
            print("[ftlk] OB2.4 probe stdout:\n", probe.stdout or "(empty)")
            raise RuntimeError("Open Babel 2.4 is not correctly configured for Py2 subprocesses")
        trg_base = os.path.join(tmpdir, 'target_cov')
        ch, rn, rs = find_residue_info(rec_pdbqt, newA[0]); cov_res = '%s:%s%s' % (ch, rn, rs)

        agfr_cmd = [
            ADFR_PY, os.path.join(PKG_BIN, 'runAGFR.py'),
            '-r', rec_pdbqt,
            '-b', 'user', '%.3f'%center[0], '%.3f'%center[1], '%.3f'%center[2],
            '%.3f'%box[0], '%.3f'%box[1], '%.3f'%box[2],
            '-o', trg_base,
            '-c', str(newA[0]), str(newA[1]),
            '-t', str(newA[2]),
            '-x', cov_res
        ]
        subprocess.check_call(agfr_cmd, env=env_py2)
        trg_file = trg_base + '.trg'
        print('[ftlk] Target created:', trg_file)

        # 6) Final docking (Py2) with the same OB2.4.1 env
        lig_for_adfr = lig_pdbqt
        base_name = os.path.basename(lig_for_adfr).replace('.pdbqt','')
        if args.jobName:
            renamed = os.path.join(tmpdir, args.jobName + '.pdbqt')
            os.rename(lig_for_adfr, renamed); lig_for_adfr = renamed; base_name = args.jobName
            print('[ftlk] Renamed temp ligand to', os.path.basename(lig_for_adfr))

        adfr_cmd = [
            ADFR_PY, os.path.join(PKG_BIN, 'runADFR.py'),
            '-l', lig_for_adfr,
            '-t', trg_file,
            '-C'
        ] + [str(x) for x in newC] + [
            '--nbRuns', str(args.nbRuns),
            '--maxEvals', str(args.maxEvals),
            '--seed', str(args.seed),
            '-O'
        ] + list(extra)
        print('[ftlk] >>', ' '.join(adfr_cmd))
        subprocess.check_call(adfr_cmd, cwd=out_dir, env=env_py2)   #

        print('[ftlk] Docking finished.')
        search_dir = out_dir if out_dir else '.'
        pat = os.path.join(search_dir, base_name + "*.dro")
        found = glob.glob(pat)
        if found:
            pattern = re.compile(rf"^{re.escape(base_name)}_.*\.dro$")
            files = os.listdir(search_dir)
            matches = [f for f in files if pattern.match(f)]
            for i in matches:
                i1=i.replace('.dro','')
                dro_path=os.path.join(search_dir,i)
                extract_path=os.path.join(search_dir,i1)
                extract_tar(dro_path,extract_path)
                path=os.path.join(search_dir, i1, f'{i1}_dro')
                ligpdbqtpath=os.path.join(path, f'{i1}_out.pdbqt')
                receptorpath=os.path.join(path, 'input', 'complex_asReceptor.pdbqt')
                #print(ligpdbqtpath,receptorpath)
                merged_path=os.path.join(path, 'merged.pdbqt')
                outpdb=os.path.join(search_dir,f'{i1}_modified.pdb')
                combine_pdbqt_convert_pdb(ligpdbqtpath, receptorpath, merged_path, outpdb)
                shutil.rmtree(extract_path)
       
            print('[ftlk] Docking result:', outpdb)
                  
        else:
            print('[ftlk] WARNING: No output matched', pat)

    except Exception:
        print('\n[ftlk] ERROR: Docking aborted:')
        traceback.print_exc()
        if not args.keepTmp:
            args.keepTmp = True
            print('[ftlk] Temp directory preserved for debugging.')
        print('[ftlk] Temp dir:', tmpdir)
        sys.exit(1)
    finally:
        if not args.keepTmp and 'tmpdir' in locals() and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == '__main__':
    main()

