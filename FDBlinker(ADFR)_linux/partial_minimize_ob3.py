#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Constrained local optimization for only selected atoms (Open Babel 3.x)
CLEAN VERSION - Relies on a correctly configured Conda environment.
"""
import os, sys
from pathlib import Path
import argparse

# --- All environment sanitation has been REMOVED. We now trust the conda env. ---
from openbabel import openbabel as ob
from openbabel import pybel

def _read_with_conv(path, fmt):
    """Helper to read a molecule using OBConversion."""
    obmol = ob.OBMol()
    conv = ob.OBConversion()
    if not conv.SetInFormat(fmt):
        return None
    if not conv.ReadFile(obmol, path) or obmol.NumAtoms() == 0:
        return None
    return obmol

def _read_obmol(path):
    """Reads a molecule file, trying common formats."""
    ext = Path(path).suffix.lstrip('.').lower()
    # Prioritize the file's extension
    trials = [ext] if ext in ("sdf", "mol2", "pdb", "mol") else []
    # Add fallbacks
    for t in ("sdf", "mol2", "pdb", "mol"):
        if t not in trials:
            trials.append(t)
    
    for t in trials:
        obmol = _read_with_conv(path, t)
        if obmol:
            return obmol
    
    # If all direct attempts fail, raise a clear error.
    raise IOError(f"Cannot read input file: {path}; tried formats: {trials}")

def _write_pdb(obmol, out_path):
    """Helper to write a molecule to a PDB file."""
    conv = ob.OBConversion()
    if not conv.SetOutFormat("pdb"):
        raise RuntimeError("Could not set PDB as output format.")
    conv.WriteFile(obmol, out_path)
    conv.CloseOutFile()

def _gen3d_preserve_fixed(obmol, fixed_idx):
    """Generates 3D coords for the whole molecule, then restores fixed atoms."""
    saved_coords = {idx: (obmol.GetAtom(idx).GetX(),
                          obmol.GetAtom(idx).GetY(),
                          obmol.GetAtom(idx).GetZ()) for idx in fixed_idx}
    
    builder = ob.OBBuilder()
    builder.Build(obmol)
    
    for idx, (x, y, z) in saved_coords.items():
        obmol.GetAtom(idx).SetVector(x, y, z)

def constrained_minimize(inp, outp, movable_indices, steps=250, ff_name='mmff94', blank_tail=0):
    """Main function for constrained minimization."""
    # Pre-read atom names from the MOL2 file to ensure they are preserved later.
    atom_names_from_mol2 = []
    try:
        with open(inp, 'r') as f:
            in_atom_block = False
            for line in f:
                if line.strip().startswith('@<TRIPOS>ATOM'):
                    in_atom_block = True
                    continue
                if line.strip().startswith('@<TRIPOS>'):
                    in_atom_block = False
                if in_atom_block and line.strip():
                    parts = line.split()
                    if len(parts) > 1:
                        atom_names_from_mol2.append(parts[1])
    except Exception as e:
        print(f"Warning: Could not pre-read atom names from {inp}. Reason: {e}", file=sys.stderr)

    obmol = _read_obmol(inp)
    natoms = obmol.NumAtoms()
    movable = set(movable_indices)
    fixed = [i for i in range(1, natoms + 1) if i not in movable]

    if any(coord.GetZ() == 0.0 for coord in ob.OBMolAtomIter(obmol)):
         _gen3d_preserve_fixed(obmol, fixed)

    ff = ob.OBForceField.FindForceField(ff_name)
    if not ff or not ff.Setup(obmol):
        print(f"Warning: Force field '{ff_name}' setup failed. Falling back to 'uff'.", file=sys.stderr)
        ff = ob.OBForceField.FindForceField('uff')
        if not ff or not ff.Setup(obmol):
            raise RuntimeError("Force field setup failed completely.")

    cons = ob.OBFFConstraints()
    for idx in fixed:
        cons.AddAtomConstraint(idx)
    ff.SetConstraints(cons)

    ff.ConjugateGradients(int(steps))
    ff.GetCoordinates(obmol)

    # Before writing to PDB, explicitly set the atom ID for each atom using
    # the names we pre-read from the original MOL2 file. This forces the
    # PDB writer to use the correct names instead of defaulting to element symbols.
    if atom_names_from_mol2 and len(atom_names_from_mol2) == obmol.NumAtoms():
        print("Restoring original atom names for PDB output...")
        for i, atom in enumerate(ob.OBMolAtomIter(obmol)):
            res = atom.GetResidue()
            if not res:
                # If an atom has no residue, create one for it (safety measure)
                res = obmol.NewResidue()
                res.AddAtom(atom)
                res.SetName("LIG")
                res.SetNum(1)
                res.SetChain("A")
            
            # The key step: Set the PDB atom ID using the preserved name
            atom_name = atom_names_from_mol2[i]
            res.SetAtomID(atom, atom_name)
    else:
        print(f"Warning: Atom count mismatch or names not read. Skipping atom name restoration.", file=sys.stderr)

    # Standardize residue info for LigandB ONLY, preserving it for the sidechain.
    # The number of sidechain atoms is passed via the '--num-new' argument (blank_tail).
    if blank_tail > 0:
        num_total_atoms = obmol.NumAtoms()
        num_ligand_b_atoms = num_total_atoms - blank_tail

        # Iterate through all atoms, but only modify the first 'num_ligand_b_atoms'
        for i, atom in enumerate(ob.OBMolAtomIter(obmol)):
            if i >= num_ligand_b_atoms:
                break  # Stop modification once we reach the sidechain atoms
            
            res = atom.GetResidue()
            if res:
                res.SetName("LIG")
                res.SetNum(1)
                res.SetChain("A")

    _write_pdb(obmol, outp)

    # Blank the names of the last N atoms if requested
    if blank_tail > 0:
        with open(outp, 'r') as f:
            lines = f.readlines()
        atom_indices = [i for i, line in enumerate(lines) if line.startswith(('ATOM', 'HETATM'))]
        for i in atom_indices[-blank_tail:]:
            line = lines[i]
            lines[i] = line[:12] + '    ' + line[16:]
        with open(outp, 'w') as f:
            f.writelines(lines)

def _parse_indices_csv(s):
    """Parses a comma-separated string of numbers into a list of ints."""
    return [int(x) for x in s.split(',') if x.strip()]

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="Constrained local optimization (Open Babel 3.x)")
    ap.add_argument('--in', dest='inp', required=True, help='Input MOL2/SDF/PDB with connectivity')
    ap.add_argument('--out', dest='outp', required=True, help='Output PDB path')
    ap.add_argument('--movable', required=True, help='CSV of 1-based indices allowed to move, e.g. "16,17,18"')
    ap.add_argument('--steps', type=int, default=250, help='Number of CG steps')
    ap.add_argument('--ff', default='mmff94', help='mmff94 or uff')
    ap.add_argument('--num-new', dest='num_new', type=int, default=0, help='Blank the PDB atom-names for the last N atoms')
    
    args = ap.parse_args()
    constrained_minimize(args.inp, args.outp, _parse_indices_csv(args.movable),
                         steps=args.steps, ff_name=args.ff, blank_tail=args.num_new)