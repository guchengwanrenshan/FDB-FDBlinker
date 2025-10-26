#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ftlk.py - Flexible two-ligand covalent docking wrapper (Python 2.7)

FINAL REVISION (Open Babel): This version now uses the Open Babel library to
perform a chemically-aware merging of the LigandA sidechain and LigandB.
This is the correct approach, as it operates on the molecular graph to ensure
the final PDBQT structure, including bonds, branches, and torsions, is valid.

PATH FIX 3 (FINAL): Correctly sets BABEL_DATADIR to the 'data' subdirectory,
which is the definitive fix for the "Cannot open element.txt" warnings.

REVISION (outputPath): Added a dedicated --outputPath argument to control
where all final output files are saved.
"""
from __future__ import print_function
import os
import sys
import argparse
import subprocess
import tempfile
import shutil
import traceback
import codecs
import re
import glob
import math
import itertools

# -----------------------------------------------------------------------------
# 1) Environment Setup
# -----------------------------------------------------------------------------
if 'ADFRHOME' not in os.environ:
    os.environ['ADFRHOME'] = os.path.abspath(os.path.dirname(sys.executable))
ADFRHOME = os.environ['ADFRHOME']

def _inject_paths():
    """
    [REVISED] Configures environment variables to ensure all necessary
    binaries, libraries, and data files for Open Babel are found.
    """
    ob_dir_name = 'OpenBabel-2.4.1'
    ob_full_path = os.path.join(ADFRHOME, ob_dir_name)
    if os.path.isdir(ob_full_path):
        print('[ftlk] Found Open Babel installation at: %s' % ob_full_path)
        
        # --- START: DEFINITIVE FIX FOR OPEN BABEL ---
        
        # 1. Set BABEL_DATADIR to the 'data' subdirectory.
        babel_data_path = os.path.join(ob_full_path, 'data')
        if os.path.isdir(babel_data_path):
            os.environ['BABEL_DATADIR'] = babel_data_path
            print('[ftlk] Set BABEL_DATADIR to: %s' % babel_data_path)
        else:
            print('[ftlk] WARNING: Open Babel "data" directory not found at: %s' % babel_data_path)
            # Fallback for older installations - set to main directory
            os.environ['BABEL_DATADIR'] = ob_full_path
            print('[ftlk] Fallback: Set BABEL_DATADIR to: %s' % ob_full_path)


        # 2. Add necessary binary/library directories to the system PATH.
        add_paths = {ob_full_path}
        for subfolder in ['bin', 'lib']:
            subfolder_path = os.path.join(ob_full_path, subfolder)
            if os.path.isdir(subfolder_path):
                add_paths.add(subfolder_path)

        cur_path_list = os.environ.get('PATH', '').split(os.pathsep)
        final_path_list = []
        # Prepend Open Babel paths to give them priority
        for p in sorted(list(add_paths)):
            if p and os.path.isdir(p) and p not in final_path_list:
                final_path_list.insert(0, p) 

        for pth in cur_path_list:
            if pth and pth not in final_path_list:
                final_path_list.append(pth)

        os.environ['PATH'] = os.pathsep.join(final_path_list)
        print('[ftlk] Environment PATH configured for Open Babel.')
        
    else:
        print('[ftlk] WARNING: Could not find the expected Open Babel directory: %s' % ob_dir_name)
    # --- END: DEFINITIVE FIX FOR OPEN BABEL ---

# Run the path injection immediately
_inject_paths()

# Now, with the environment fixed, we can attempt the import.
try:
    import pybel
except ImportError:
    print("\n[ftlk] FATAL ERROR: The 'openbabel' Python library could not be imported.")
    print("[ftlk] Please check your ADFRsuite installation and PATH configuration.")
    sys.exit(1)


if sys.platform.startswith('win'):
    sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
    sys.stderr = codecs.getwriter('utf-8')(sys.stderr)

# -----------------------------------------------------------------------------
# The rest of the script remains unchanged.
# -----------------------------------------------------------------------------
PYTHON_EXE = 'python'
PKG_BIN = os.path.join(ADFRHOME, 'Lib', 'site-packages', 'ADFR', 'bin')
_RX_SKIP = re.compile(r'^(ROOT|ENDROOT|BRANCH|ENDBRANCH|TORSDOF|REMARK)')

def create_complex_receptor_file(original_pdb, final_receptor_path):
    print('[ftlk] Preparing complex file...')
    workdir = os.path.dirname(final_receptor_path)
    protein_pdb_path = os.path.join(workdir, 'protein_original.pdb')
    ligand_pdb_path = os.path.join(workdir, 'ligandA_original.pdb')

    with open(original_pdb, 'r') as fin, \
         open(protein_pdb_path, 'w') as fp, \
         open(ligand_pdb_path, 'w') as fl:
        for line in fin:
            if line.startswith('ATOM'):
                fp.write(line)
            elif line.startswith('HETATM'):
                fl.write(line)
            elif line.startswith('TER'):
                fp.write(line)

    if not os.path.getsize(ligand_pdb_path) > 0:
        raise RuntimeError('HETATM records for LigandA not found in %s.' % (original_pdb))

    mgl_util = os.path.join(ADFRHOME, 'Lib', 'site-packages', 'AutoDockTools', 'Utilities24')
    prep_rec_script = os.path.join(mgl_util, 'prepare_receptor4.py')
    prep_lig_script = os.path.join(mgl_util, 'prepare_ligand4.py')
    protein_qt_path = protein_pdb_path.replace('.pdb', '.pdbqt')
    ligand_qt_path = ligand_pdb_path.replace('.pdb', '.pdbqt')

    subprocess.check_call([PYTHON_EXE, prep_rec_script, '-r', os.path.basename(protein_pdb_path), '-A', 'checkhydrogens', '-U', 'nphs_lps_waters', '-o', os.path.basename(protein_qt_path)], cwd=workdir)
    subprocess.check_call([PYTHON_EXE, prep_lig_script, '-l', os.path.basename(ligand_pdb_path), '-A', 'checkhydrogens', '-o', os.path.basename(ligand_qt_path)], cwd=workdir)

    final_lines, total_atom_count = [], 0
    last_protein_atom_line = None

    with open(protein_qt_path, 'r') as f_prot:
        for line in f_prot:
            if _RX_SKIP.match(line) or line.startswith('TER'):
                continue
            if line.startswith('ATOM'):
                total_atom_count += 1
                final_lines.append('%s%5d%s' % (line[:6], total_atom_count, line[11:]))
                last_protein_atom_line = line

    if last_protein_atom_line:
        res_name = last_protein_atom_line[17:20]
        chain_id = last_protein_atom_line[21:22]
        res_seq = last_protein_atom_line[22:26]
        ter_serial = total_atom_count + 1
        final_lines.append('TER   %5d      %s %s%s\n' % (ter_serial, res_name, chain_id, res_seq))
        total_atom_count += 1

    with open(ligand_qt_path, 'r') as f_lig:
        for line in f_lig:
            if _RX_SKIP.match(line):
                continue
            if line.startswith('HETATM'):
                total_atom_count += 1
                final_lines.append('ATOM  %5d%s' % (total_atom_count, line[11:]))

    with open(final_receptor_path, 'w') as fout:
        fout.write("".join(final_lines))
    print('[ftlk] Final receptor file written to %s (%d atoms)' % (final_receptor_path, total_atom_count))

def find_new_anchor_serials(original_pdb_path, new_pdbqt_path, original_serial_list):
    new_serials = []

    with open(new_pdbqt_path, 'r') as f_new:
        lines = f_new.readlines()

    for orig_serial_to_find in original_serial_list:
        atom_name_to_find, resname_to_find, resseq_to_find = None, None, None
        with open(original_pdb_path, 'r') as f_orig:
            for line in f_orig:
                if line.startswith('HETATM'):
                    try:
                        if int(line[6:11]) == orig_serial_to_find:
                            atom_name_to_find = line[12:16].strip()
                            resname_to_find = line[17:20].strip()
                            resseq_to_find = line[22:26].strip()
                            break
                    except (ValueError, IndexError):
                        continue

        if atom_name_to_find is None:
            raise RuntimeError("Could not find HETATM record for original serial %d in %s." % (orig_serial_to_find, original_pdb_path))

        found_new_serial = None
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    if (line[12:16].strip() == atom_name_to_find and
                        line[17:20].strip() == resname_to_find and
                        line[22:26].strip() == resseq_to_find):
                        found_new_serial = int(line[6:11])
                        break
                except (ValueError, IndexError):
                    continue

        if found_new_serial is None:
            raise RuntimeError("Could not find atom '%s' in residue '%s%s' in '%s'." % (atom_name_to_find, resname_to_find, resseq_to_find, new_pdbqt_path))

        new_serials.append(found_new_serial)

    return new_serials

def get_atom_coords_by_serial(pdbqt_path, serial_number):
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and int(line[6:11]) == serial_number:
                return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    raise ValueError('Atom with serial number %d not found in %s.' % (serial_number, pdbqt_path))

def find_residue_info(pdbqt_path, atom_serial):
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and int(line[6:11]) == atom_serial:
                chain_id = line[21:22].strip()
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                if not chain_id: chain_id = ' '
                return chain_id, res_name, res_seq
    raise ValueError('Atom with serial number %d not found in %s.' % (atom_serial, pdbqt_path))

def _get_atom_info_manually(pdbqt_path, serial_to_find):
    """
    Manually parses a PDBQT file to get reliable atom info,
    bypassing pybel's faulty PDBQT reader.
    """
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    current_serial = int(line[6:11])
                    if current_serial == serial_to_find:
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        chain_id = line[21:22].strip()
                        res_num_str = line[22:26].strip()
                        
                        return {
                            'name': atom_name,
                            'res_name': res_name,
                            'chain': chain_id,
                            'res_num': int(res_num_str)
                        }
                except (ValueError, IndexError):
                    continue
    raise ValueError("Atom with serial %d not found in %s during manual parsing." % (serial_to_find, pdbqt_path))

def _correct_residue_info_in_pdbqt(pdbqt_path, sidechain_atom_infos):
    """
    Manually corrects the residue information in the final PDBQT file.
    This is a post-processing step after prepare_ligand4.py runs.
    """
    print('[ftlk] Post-processing: Restoring original residue info in %s...' % os.path.basename(pdbqt_path))
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()

    atom_lines_indices = [i for i, line in enumerate(lines) if line.startswith(('ATOM', 'HETATM'))]
    
    num_sidechain_atoms = len(sidechain_atom_infos)
    if len(atom_lines_indices) < num_sidechain_atoms:
        print("[ftlk] WARNING: Not enough atoms in PDBQT to restore residue info. Skipping.")
        return

    # Target the last N atoms, which are assumed to be the sidechain
    target_indices = atom_lines_indices[-num_sidechain_atoms:]

    if len(target_indices) != len(sidechain_atom_infos):
         print("[ftlk] WARNING: Mismatch between target atoms and info records. Aborting residue correction.")
         return

    for i in range(num_sidechain_atoms):
        line_idx = target_indices[i]
        info = sidechain_atom_infos[i]
        line = lines[line_idx]

        # PDB format for residue info: resname (3), chain (1), resnum (4)
        # Example: 'CYS A 164'
        res_name = info['res_name'][:3].strip()
        chain_id = info['chain'][:1].strip()
        res_num_str = str(info['res_num'])

        # Format the residue segment: "CYS A 164 "
        # line[17:27] is the slice for "RES C RESI "
        res_segment = "{:>3s} {:1s}{:>4s} ".format(res_name, chain_id, res_num_str)
        
        # Replace the segment in the line
        lines[line_idx] = line[:17] + res_segment + line[27:]

    with open(pdbqt_path, 'w') as f:
        f.writelines(lines)

def _fix_sidechain_info_in_pdb(file_path, sidechain_atom_names, sidechain_atom_infos):
    """
    Corrects sidechain info based on blank names, and then
    forcefully converts ALL 'ATOM' records in the file to 'HETATM' to ensure
    maximum compatibility with prepare_ligand4.py.
    """
    print('[ftlk] Polishing file: Correcting sidechain atom names and residue info in %s...' % os.path.basename(file_path))
    with open(file_path, 'r') as f:
        lines = f.readlines()

    num_sidechain_atoms = len(sidechain_atom_names)

    # --- Part 1: Find and correct the specific sidechain atoms by blank name ---
    target_line_indices = []
    all_atom_line_indices = [i for i, line in enumerate(lines) if line.startswith(('ATOM', 'HETATM'))]

    for line_idx in all_atom_line_indices:
        line = lines[line_idx]
        if not line[12:16].strip():
            target_line_indices.append(line_idx)
    
    if len(target_line_indices) != num_sidechain_atoms:
        print("[ftlk] ERROR: Expected to find %d atoms with blank names, but found %d. Cannot proceed with corrections." % (num_sidechain_atoms, len(target_line_indices)))
        raise RuntimeError("Sidechain correction failed due to mismatch in expected blank-named atoms.")

    print('[ftlk] Found %d target atoms. Applying specific corrections...' % len(target_line_indices))

    for i in range(num_sidechain_atoms):
        line_idx = target_line_indices[i]
        line = lines[line_idx]

        # 1. Correct the Atom Name
        atom_name_to_fix = sidechain_atom_names[i]
        formatted_name = (" " + atom_name_to_fix).ljust(4)
        line = line[:12] + formatted_name + line[16:]

        # 2. Correct the Residue Info
        info = sidechain_atom_infos[i]
        res_name = info['res_name'][:3].strip()
        chain_id = info['chain'][:1].strip()
        res_num_str = str(info['res_num'])
        res_segment = "{:>3s} {:1s}{:>4s} ".format(res_name, chain_id, res_num_str)
        line = line[:17] + res_segment + line[27:]
        
        lines[line_idx] = line

    # --- Part 2: Force ALL ATOM records in the entire file to HETATM ---
    print('[ftlk] Forcing all ATOM records to HETATM for compatibility...')
    conversion_count = 0
    for i, line in enumerate(lines):
        if line.startswith('ATOM  '):
            lines[i] = "HETATM" + line[6:]
            conversion_count += 1
    print('[ftlk] Converted %d remaining ATOM records to HETATM.' % conversion_count)

    with open(file_path, 'w') as f:
        f.writelines(lines)

def _restore_atom_record_type_in_pdbqt(pdbqt_path, num_atoms_to_change=3):
    """
    Post-processes the final PDBQT file to change the record type of the
    first N atoms from HETATM back to ATOM. This restores the chemical
    distinction of the sidechain atoms after prepare_ligand4.py has run.
    """
    print('[ftlk] Post-processing: Restoring ATOM records for sidechain in %s...' % os.path.basename(pdbqt_path))
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()

    atom_lines_found = 0
    for i, line in enumerate(lines):
        if atom_lines_found >= num_atoms_to_change:
            break # Stop after we've processed the required number of atoms

        if line.startswith('HETATM'):
            lines[i] = "ATOM  " + line[6:]
            atom_lines_found += 1
    
    if atom_lines_found < num_atoms_to_change:
        print("[ftlk] WARNING: Expected to change %d records to ATOM, but only found %d HETATM records at the start of the file." % (num_atoms_to_change, atom_lines_found))

    with open(pdbqt_path, 'w') as f:
        f.writelines(lines)

def _swap_root_atoms_and_fix_branch(pdbqt_path):
    """
    Swaps the first two atoms in the ROOT block and corrects
    the subsequent branch definition, as per user request.
    This is a specific post-processing step to match a desired output format.
    Specifically:
    1. Swaps the line content of ATOM 1 and ATOM 2 (from column 12 onwards).
    2. Finds the first "BRANCH   1  ..." line and changes it to "BRANCH   2  ...".
    """
    print('[ftlk] Applying user-requested post-processing: Swapping root atoms and fixing branch.')
    
    try:
        with open(pdbqt_path, 'r') as f:
            lines = f.readlines()

        # Find the line index of the 'ROOT' keyword
        root_idx = -1
        for i, line in enumerate(lines):
            if line.startswith('ROOT'):
                root_idx = i
                break
        
        if root_idx == -1:
            print("[ftlk] WARNING: 'ROOT' keyword not found. Skipping atom swap fix.")
            return

        atom1_line_idx = root_idx + 1
        atom2_line_idx = root_idx + 2

        line1 = lines[atom1_line_idx]
        line2 = lines[atom2_line_idx]

        # Ensure these two lines are ATOM/HETATM records
        if (line1.startswith(('ATOM', 'HETATM'))) and (line2.startswith(('ATOM', 'HETATM'))):
            # In PDB/PDBQT format, the atom serial number is in columns 7-11. 
            # We swap all content after column 11.
            content1 = line1[11:]
            content2 = line2[11:]
            
            lines[atom1_line_idx] = line1[:11] + content2
            lines[atom2_line_idx] = line2[:11] + content1
            
            print('[ftlk] -> Swapped content of ATOM 1 and ATOM 2.')

            # Find the line index of the 'ENDROOT' keyword
            endroot_idx = -1
            for i in range(root_idx, len(lines)):
                 if lines[i].startswith('ENDROOT'):
                     endroot_idx = i
                     break
            
            if endroot_idx != -1:
                # Search for the first BRANCH line after ENDROOT and correct it
                for i in range(endroot_idx + 1, len(lines)):
                    line_to_check = lines[i]
                    if line_to_check.startswith('BRANCH'):
                        # The AutoDockTools PDBQT format is "BRANCH%4d%4d"
                        try:
                            first_atom_serial = int(line_to_check[6:10].strip())
                            if first_atom_serial == 1:
                                # Rebuild the line, changing the first atom serial from 1 to 2, 
                                # while keeping the rest of the line intact.
                                rest_of_line = line_to_check[10:]
                                lines[i] = "BRANCH%4d%s" % (2, rest_of_line)
                                print('[ftlk] -> Corrected BRANCH definition from atom 1 to 2.')
                                break # Only correct the first one found
                        except (ValueError, IndexError):
                            # If the line format is incorrect, skip it
                            continue
        else:
            print("[ftlk] WARNING: Expected ATOM/HETATM lines after 'ROOT' not found. Skipping atom swap.")

        with open(pdbqt_path, 'w') as f:
            f.writelines(lines)

    except (IndexError, ValueError) as e:
        print("[ftlk] WARNING: Could not perform root atom swap due to unexpected file format. Error: %s" % e)

#
# NEW HELPER FUNCTION TO PRE-REPAIR A PDB FILE (Definitive Version with Manual PDB Writer)
#
def _pre_repair_pdb(input_pdb_path, output_pdb_path):
    """
    Uses AD4LigandPreparation to perform a standardized repair on a PDB file.
    This function adds bonds, hydrogens, and Gasteiger charges, then saves
    the result as a new, chemically complete PDB file. This stabilizes the
    molecule for all subsequent processing steps.
    
    Args:
        input_pdb_path (str): Path to the input PDB file.
        output_pdb_path (str): Path where the repaired PDB file will be saved.
    """
    # Import necessary classes from MGLTools.
    from MolKit import Read
    from AutoDockTools.MoleculePreparation import AD4LigandPreparation
    # --- THIS IS THE CORRECTED LINE ---
    from MolKit.pdbWriter import PdbWriter

    print('[ftlk_repair] Reading molecule from: %s' % os.path.basename(input_pdb_path))
    mols = Read(input_pdb_path)
    if not mols:
        raise ValueError("MGLTools could not read molecule from %s" % input_pdb_path)
    
    mol = mols[0]
    if len(mols) > 1:
        print('[ftlk_repair] Multiple molecules found; selecting the one with the most atoms.')
        for m in mols[1:]:
            if len(m.allAtoms) > len(mol.allAtoms):
                mol = m

    print('[ftlk_repair] Performing repairs: building bonds, adding hydrogens, and calculating charges...')
    
    preparer = AD4LigandPreparation(
        mol,
        mode='automatic',
        repairs='bonds_hydrogens',
        charges_to_add='gasteiger',
        cleanup='',
        allowed_bonds="backbone",
        root='auto',
        outputfilename=None
    )

    print('[ftlk_repair] Writing repaired molecule to: %s' % os.path.basename(output_pdb_path))
    writer = PdbWriter()
    writer.write(output_pdb_path, mol, records=['ATOM', 'HETATM', 'CONECT'])



def create_merged_ligand_pdbqt(ligand_b_path, complex_path, atoms_a_serials, atom_b_serial, final_output_path, workdir, python_exe):
    """
    Creates a valid PDBQT file for covalent docking.
    This version pre-simulates prepare_ligand4.py's cleanup process to find the
    correct, post-cleanup root index, resolving the IndexError.
    """
    print('[ftlk] Starting robust ligand creation workflow with root selection...')
    
    # --- Part 1 & 2: Build and Unify Molecule in Memory ---
    print('[ftlk] [Step 1/5] Building merged molecule and collecting original residue data...')
    try:
        ligand_b_mol = next(pybel.readfile("pdbqt", ligand_b_path))
    except StopIteration:
        raise IOError("Open Babel could not read LigandB PDBQT file.")

    sidechain_mol = pybel.ob.OBMol()
    element_table = pybel.ob.OBElementTable()
    map_serial_to_new_idx = {}
    sidechain_names_in_order = []
    sidechain_original_infos = []

    for serial in atoms_a_serials:
        atom_info = _get_atom_info_manually(complex_path, serial)
        sidechain_original_infos.append(atom_info)
        atom_name = atom_info['name']
        sidechain_names_in_order.append(atom_name) 
        
        element_symbol = atom_name[0].upper()
        atomic_num = element_table.GetAtomicNum(element_symbol)
        new_atom = sidechain_mol.NewAtom()
        new_atom.SetAtomicNum(atomic_num)
        
        new_residue = sidechain_mol.NewResidue()
        new_residue.SetName(atom_info['res_name'])
        new_residue.SetNum(atom_info['res_num'])
        new_residue.SetChain(atom_info['chain'])
        new_residue.AddAtom(new_atom)
        # FIX: Explicitly set the PDB-style atom name for the new atom.
        # The .encode('utf-8') is for compatibility with the wrapper.
        new_residue.SetAtomID(new_atom, atom_name.encode('utf-8'))

        map_serial_to_new_idx[serial] = new_atom.GetIdx()
    
    if len(atoms_a_serials) > 1:
        sidechain_mol.AddBond(map_serial_to_new_idx[atoms_a_serials[1]], map_serial_to_new_idx[atoms_a_serials[0]], 1)
    if len(atoms_a_serials) > 2:
        sidechain_mol.AddBond(map_serial_to_new_idx[atoms_a_serials[0]], map_serial_to_new_idx[atoms_a_serials[2]], 1)

    num_atoms_in_ligand_b = ligand_b_mol.OBMol.NumAtoms()
    ligand_b_mol.OBMol += sidechain_mol
    
    attach_b_idx = ligand_b_mol.atoms[atom_b_serial - 1].idx
    attach_a_idx = map_serial_to_new_idx[atoms_a_serials[1]] + num_atoms_in_ligand_b
    ligand_b_mol.OBMol.AddBond(attach_b_idx, attach_a_idx, 1)

    print('[ftlk] [Step 2/5] Temporarily unifying all residue info to "LIG" for compatibility...')
    for atom in pybel.ob.OBMolAtomIter(ligand_b_mol.OBMol):
        res = atom.GetResidue()
        if res:
            res.SetName("LIG")
            res.SetNum(1)
            res.SetChain("A")

    # --- Part 3 (REPLACED): Constrained 3D build & local optimization using Python 3.8 + OB3 ---
    intermediate_mol2_path = os.path.join(workdir, 'merged_ligand_tmp.mol2')
    ligand_b_mol.write("mol2", intermediate_mol2_path, overwrite=True)

    # 2) Call external Python 3.8 helper using Open Babel 3.0.0 to ONLY optimize the 3 new atoms
    python38 = r"C:\Users\ROG\anaconda3\envs\falk_env\python.exe"  # your Python 3.8
    helper_py = r"C:\Program Files (x86)\ADFRsuite-1.0\Lib\site-packages\ADFR\bin\partial_minimize_ob3.py" # path to the new helper

    # The sidechain atoms were appended after LigandB atoms:
    movable_indices = [num_atoms_in_ligand_b + i for i in range(1, len(atoms_a_serials)+1)]
    movable_csv = ",".join(str(i) for i in movable_indices)

    intermediate_pdb_path = os.path.join(workdir, 'merged_ligand.pdb')  # downstream expects PDB

    clean_env = os.environ.copy()

    orig_path = clean_env.get('PATH', '')
    kept = [p for p in orig_path.split(os.pathsep) if "OpenBabel-2.4.1" not in p]

    env_root = os.path.dirname(python38)
    env_scripts = os.path.join(env_root, 'Scripts')
    env_libbin = os.path.join(env_root, 'Library', 'bin')
    prefixes = [p for p in (env_libbin, env_scripts, env_root) if os.path.isdir(p)]

    clean_env['PATH'] = os.pathsep.join(prefixes + kept)

    for cand in (os.path.join(env_root, 'share', 'openbabel'),
                os.path.join(env_root, 'Library', 'share', 'openbabel')):
        if os.path.isdir(cand):
            clean_env['BABEL_DATADIR'] = cand
            break

    for cand in (os.path.join(env_root, 'Library', 'lib', 'openbabel', '3.1'),
                os.path.join(env_root, 'Library', 'lib', 'openbabel', '3.1.1'),
                os.path.join(env_root, 'lib', 'openbabel', '3.1'),
                os.path.join(env_root, 'lib', 'openbabel', '3.1.1')):
        if os.path.isdir(cand):
            clean_env['OB_PLUGIN_PATH'] = cand
            break

    for k in list(clean_env.keys()):
        if k.upper().startswith('BABEL_') and k != 'BABEL_DATADIR':
            del clean_env[k]


    # Keep the last N atom names blank so _fix_sidechain_info_in_pdb() works unchanged
    subprocess.check_call([
        python38, helper_py,
        "--in", intermediate_mol2_path,
        "--out", intermediate_pdb_path,
        "--movable", movable_csv,
        "--steps", "250",
        "--ff", "mmff94",
        "--num-new", str(len(atoms_a_serials))
    ], env=clean_env) # <-- Pass the cleaned environment to the subprocess


    _fix_sidechain_info_in_pdb(
        intermediate_pdb_path, 
        sidechain_names_in_order, 
        sidechain_original_infos
    )
    
    repaired_pdb_path = os.path.join(workdir, 'merged_ligand_repaired.pdb')
    
    print('\n[ftlk] STEP 3.5: Pre-repairing the merged PDB to ensure chemical completeness.')
    _pre_repair_pdb(intermediate_pdb_path, repaired_pdb_path)
    print('[ftlk] Pre-repair complete. All subsequent steps will use the repaired file.\n')

    print('[ftlk] Pre-simulating cleanup to find the correct post-cleanup root index...')
    from MolKit import Read
    from AutoDockTools.MoleculePreparation import AutoDockMoleculePreparation
    
    # MolKit's Read needs to be in the correct directory
    mol = Read(repaired_pdb_path)[0] 
    mol.buildBondsByDistance()

    # Perform the exact same cleanup that prepare_ligand4.py does by default.
    # We set charges_to_add=None because we only care about the cleanup process 
    # which affects atom count and order, not the charges at this stage.
    AutoDockMoleculePreparation(
        mol,
        mode='automatic',
        repairs='',
        charges_to_add='gasteiger',
        cleanup='nphs_lps', # Default cleanup for ligands
        outputfilename=None,
        debug=False,
        version=4
    )

    # In the cleaned mol.allAtoms, find our target root atom.
    # The original serial was num_atoms_in_ligand_b + len(atoms_a_serials) -> 15 + 3 = 18.
    target_serial_pre_cleanup = num_atoms_in_ligand_b + len(atoms_a_serials)
    root_index_post_cleanup = None
    
    # Find the atom that corresponds to our original TORSION_ANCHOR.
    # We use its name ('C2') and original residue name ('l01') as identifiers.
    torsion_anchor_info = sidechain_original_infos[2] # 0:BOND_ATOM_1, 1:BOND_ATOM_2, 2:TORSION_ANCHOR
    torsion_anchor_name = sidechain_names_in_order[2]

    for i, at in enumerate(mol.allAtoms):
        # We need to check against the corrected residue info
        resname_in_mol = at.parent.type if hasattr(at, 'parent') else ''
        if at.name == torsion_anchor_name and resname_in_mol == torsion_anchor_info['res_name']:
            root_index_post_cleanup = i
            break

    if root_index_post_cleanup is None:
        raise RuntimeError("Could not locate the root atom (%s in %s) after cleanup simulation." % (torsion_anchor_name, torsion_anchor_info['res_name']))

    print('[ftlk] Determined pre-cleanup root serial was %d. Post-cleanup root index is %d.' % (target_serial_pre_cleanup, root_index_post_cleanup))


    # --- Part 4: Generate PDBQT using the correct, post-cleanup root index ---
    print('[ftlk] [Step 4/5] Creating PDBQT with prepare_ligand4.py, using correct root index...')
    mgl_util = os.path.join(ADFRHOME, 'Lib', 'site-packages', 'AutoDockTools', 'Utilities24')
    prep_lig_script = os.path.join(mgl_util, 'prepare_ligand4.py')
    
    cmd_prepare_ligand = [
        python_exe, prep_lig_script,
        '-l', os.path.basename(repaired_pdb_path),
        '-o', os.path.basename(final_output_path),
        '-A', 'checkhydrogens',
        '-R', str(root_index_post_cleanup), # Pass the correct, post-cleanup index
        '-v' # Add verbose flag for better debugging if issues persist
    ]
    subprocess.check_call(cmd_prepare_ligand, cwd=workdir)
    print('[ftlk] Intermediate PDBQT created at: %s' % final_output_path)
    
    # --- Part 5: POST-PROCESS the final PDBQT to restore ATOM records ---
    print('[ftlk] [Step 5/5] Finalizing file: Restoring ATOM records for sidechain atoms...')
    _restore_atom_record_type_in_pdbqt(final_output_path, num_atoms_to_change=len(atoms_a_serials))
    _swap_root_atoms_and_fix_branch(final_output_path)

    new_covalent_serials = [1, 2, 3] 
    
    return new_covalent_serials

def process_model(model_lines):
    """
    Processes a list of lines corresponding to a single MODEL...ENDMDL block.
    This function is called by run_postprocessing.
    """
    processed_lines = []
    atom_map = {}  # Maps old atom serials to new ones

    # Step 1: Find all atoms to keep (serial > 3) and build the serial map
    atoms_to_keep = []
    for line in model_lines:
        if line.startswith(('ATOM', 'HETATM')):
            try:
                serial = int(line[6:11])
                if serial > 3:
                    atoms_to_keep.append(line)
            except (ValueError, IndexError):
                continue

    if not atoms_to_keep:
        for line in model_lines:
            if not line.startswith(('ATOM', 'HETATM', 'ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
                processed_lines.append(line)
        return processed_lines

    for i, line in enumerate(atoms_to_keep):
        original_serial = int(line[6:11])
        atom_map[original_serial] = i + 1

    # Step 2: Iterate through original lines and build the new, valid structure
    is_root_placed = False
    for line in model_lines:
        if line.startswith('MODEL') or line.startswith('USER') or line.startswith('REMARK'):
            processed_lines.append(line)
            continue

        if line.startswith(('ATOM', 'HETATM')):
            try:
                old_serial = int(line[6:11])
                if old_serial in atom_map:
                    new_serial = atom_map[old_serial]
                    new_line = "{:6s}{:5d}{}".format(line[:6], new_serial, line[11:])
                    if not is_root_placed:
                        processed_lines.append("ROOT\n")
                        processed_lines.append(new_line)
                        processed_lines.append("ENDROOT\n")
                        is_root_placed = True
                    else:
                        processed_lines.append(new_line)
            except (ValueError, IndexError):
                continue

        elif line.startswith(('BRANCH', 'ENDBRANCH')):
            try:
                parts = line.split()
                keyword = parts[0]
                old_s1, old_s2 = int(parts[1]), int(parts[2])
                if old_s1 in atom_map and old_s2 in atom_map:
                    new_s1, new_s2 = atom_map[old_s1], atom_map[old_s2]
                    processed_lines.append("{:6s}{:4d}{:4d}\n".format(keyword, new_s1, new_s2))
            except (ValueError, IndexError):
                continue
        
        elif line.startswith('ENDMDL'):
            processed_lines.append(line)

    return processed_lines

def run_postprocessing(input_file, output_file):
    """
    Reads a PDBQT, processes each model to remove the sidechain, and writes a new file.
    """
    print('[ftlk] Starting post-processing of %s...' % input_file)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        model_chunk = []
        for line in infile:
            if line.startswith('MODEL'):
                if model_chunk:
                    processed = process_model(model_chunk)
                    outfile.writelines(processed)
                model_chunk = [line]
            elif line.startswith('ENDMDL'):
                if model_chunk:
                    model_chunk.append(line)
                    processed = process_model(model_chunk)
                    outfile.writelines(processed)
                    model_chunk = []
            elif model_chunk:
                model_chunk.append(line)
            else:
                outfile.write(line)
        
        if model_chunk:
            processed = process_model(model_chunk)
            outfile.writelines(processed)
    print('[ftlk] Post-processing complete. Cleaned file saved to: %s' % output_file)

def write_combined_complex(original_complex_pdb, docked_poses_pdbqt, final_complex_pdb, sidechain_serials_to_replace, cov_chain_id, cov_res_name, cov_res_seq):
    """
    [DEFINITIVE VERSION 2.0] Combines the protein, the unmodified part of LigandA, and each docked
    ligand pose.
    This version includes final micro-adjustments for atom name and end-of-line column alignment
    to precisely match the user's reference format.
    """
    print('[ftlk] Creating final docked complex file (adjusting final column alignment)...')

    serials_to_replace_set = set(sidechain_serials_to_replace)

    static_lines = []
    with open(original_complex_pdb, 'r') as f_original:
        for line in f_original:
            if line.startswith(('ATOM', 'TER')):
                static_lines.append(line)
            elif line.startswith('HETATM'):
                try:
                    atom_serial = int(line[6:11])
                    if atom_serial not in serials_to_replace_set:
                        static_lines.append(line)
                except (ValueError, IndexError):
                    static_lines.append(line)

    with open(docked_poses_pdbqt, 'r') as f_dock, open(final_complex_pdb, 'w') as f_out:

        def process_and_write_model(ligand_atom_lines, model_tag):
            if not ligand_atom_lines:
                return

            f_out.write(model_tag)
            atom_counter = 0

            for line in static_lines:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_counter += 1
                    f_out.write('%s%5d%s' % (line[:6], atom_counter, line[11:]))
                else:
                    f_out.write(line)

            for line in ligand_atom_lines:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_counter += 1
                    
                    # --- START: FINAL ALIGNMENT FIXES ---

                    # 1. Re-format the atom name to be centered, fixing the 'S1' alignment.
                    #    This takes the atom name like 'S1', removes spaces, and centers it in a 4-char field.
                    original_atom_name = line[12:16].strip()
                    formatted_atom_name = original_atom_name.center(4) # 'S1' -> ' S1 '

                    # 2. Rebuild the line prefix with the correctly aligned atom name
                    new_res_info = "{:>3s} {:1s}{:>4s}".format(cov_res_name[:3], cov_chain_id, cov_res_seq.strip())
                    line_prefix = 'HETATM%5d %s' % (atom_counter, formatted_atom_name) # Note the space after %5d
                    
                    line_suffix = line[26:]
                    modified_line = line_prefix + new_res_info + line_suffix
                    
                    main_part = modified_line[:66]
                    
                    # 3. Derive the element symbol
                    ad_atom_type = line.strip().split()[-1]
                    if ad_atom_type == 'A':
                        element_symbol = 'C'
                    else:
                        element_symbol = ad_atom_type[0]
                    
                    # 4. The "second to last column" is the duplicated Chain ID
                    second_to_last_col = cov_chain_id
                    
                    # 5. Assemble the final line with corrected end-of-line spacing.
                    #    The change from %-70s to %-71s adds the extra space you noted.
                    final_line = '%-71s%2s   %2s\n' % (main_part, second_to_last_col.strip(), element_symbol.strip())
                    # --- END: FINAL ALIGNMENT FIXES ---
                    
                    f_out.write(final_line)
            
            f_out.write("ENDMDL\n")

        ligand_atoms_in_current_model = []
        model_tag_line = 'MODEL        1\n'

        for line in f_dock:
            if line.startswith('MODEL'):
                process_and_write_model(ligand_atoms_in_current_model, model_tag_line)
                ligand_atoms_in_current_model = []
                model_tag_line = line
            
            elif line.startswith(('ATOM', 'HETATM')):
                ligand_atoms_in_current_model.append(line)

        if ligand_atoms_in_current_model:
            process_and_write_model(ligand_atoms_in_current_model, model_tag_line)
    
    print('[ftlk] Final complex file created successfully.')

def calculate_ligand_diameter(ligand_pdbqt_path):
    """
    Calculates the maximum distance between any two atoms in a ligand.
    This is used to determine the molecule's diameter.
    
    Args:
        ligand_pdbqt_path (str): The file path to the ligand's PDBQT file.
        
    Returns:
        float: The maximum inter-atomic distance in Angstroms.
    """
    print('[ftlk] Reading LigandB to determine its size: %s' % ligand_pdbqt_path)
    try:
        mol = next(pybel.readfile("pdbqt", ligand_pdbqt_path))
    except StopIteration:
        raise IOError("Open Babel could not read LigandB file for size calculation.")
    
    coords = [atom.coords for atom in mol.atoms]
    
    if len(coords) < 2:
        return 0.0 # Not enough atoms to calculate a distance

    max_dist_sq = 0.0
    # Use itertools.combinations to efficiently get all unique pairs of atoms
    for coord1, coord2 in itertools.combinations(coords, 2):
        dist_sq = (coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2
        if dist_sq > max_dist_sq:
            max_dist_sq = dist_sq
            
    return math.sqrt(max_dist_sq)

def main():
    parser = argparse.ArgumentParser(
        description='Two-ligand covalent docking wrapper using a two-pass AGFR strategy.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--complex', required=True, help='Protein+LigandA complex (PDB ONLY)')
    parser.add_argument('--ligand', required=True, help='LigandB to be docked (PDBQT)')
    parser.add_argument('--atomsA', nargs=3, type=int, required=True, metavar=('BOND_ATOM_1', 'BOND_ATOM_2', 'TORSION_ANCHOR'), 
                        help='Three anchor atom indices on LigandA (from original PDB).\n'
                             'Order is crucial: B-A1, A1-A2, A2-A3 bonds will be formed.')
    parser.add_argument('--atomsB', type=int, required=True, metavar=('ATTACH_ATOM'), help='Single attachment atom index on LigandB (1-based, from PDBQT file)')
    parser.add_argument('--box_size', type=float, default=None, help='Size of the cubic docking box in Angstroms. If not provided, it will be auto-calculated (LigandB diameter * sqrt(2)).')
    parser.add_argument('--nbRuns', type=int, default=8, help='Number of GA runs')
    parser.add_argument('--maxEvals', type=int, default=250000, help='Max evaluations per run')
    parser.add_argument('--seed', type=int, default=1, help='Random seed')
    parser.add_argument('--keepTmp', action='store_true', help='Keep the temp directory')
    parser.add_argument('--jobName', type=str, default=None, help='User-defined name for output files')
    parser.add_argument('--outputPath', type=str, default=None, help='Directory to save all output files.')
    args, extra = parser.parse_known_args()
    
    if not args.complex.lower().endswith('.pdb'):
        print("[ftlk] ERROR: The --complex argument must be a PDB file for this workflow.")
        sys.exit(1)

    output_dir = args.outputPath
    if output_dir:
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
                print('[ftlk] Created output directory: %s' % output_dir)
            except OSError as e:
                print('[ftlk] ERROR: Could not create output directory %s. Error: %s' % (output_dir, e))
                sys.exit(1)
        print('[ftlk] Output will be saved to: %s' % os.path.abspath(output_dir))

    tmpdir = tempfile.mkdtemp(prefix='ftlk_')
    print('[ftlk] Temp dir:', tmpdir)

    try:
        final_rec_file = os.path.join(tmpdir, 'complex_asReceptor.pdbqt')
        create_complex_receptor_file(args.complex, final_rec_file)
        
        new_atomsA_serials = find_new_anchor_serials(args.complex, final_rec_file, args.atomsA)
        print('[ftlk] Translated LigandA anchor atoms from original PDB serials %s -> PDBQT serials %s' % (args.atomsA, new_atomsA_serials))

        ligandB_with_sidechain_path = os.path.join(tmpdir, 'ligandB_with_sidechain.pdbqt')
        new_covalent_atom_serials_for_adfr = create_merged_ligand_pdbqt(
            args.ligand,
            final_rec_file,
            new_atomsA_serials,
            args.atomsB,
            ligandB_with_sidechain_path,
            tmpdir, 
            PYTHON_EXE
        )
        print('[ftlk] New covalent atom serials for ADFR -C option are: %s' % new_covalent_atom_serials_for_adfr)

        print('[ftlk] Defining docking box...')
        center = get_atom_coords_by_serial(final_rec_file, new_atomsA_serials[0])
        box_dim = 0.0
        if args.box_size is not None:
            box_dim = args.box_size
            print('[ftlk] Using user-defined box size: %.3f Angstroms' % box_dim)
        else:
            print('[ftlk] Automatically calculating box size based on LigandB dimensions...')
            ligand_diameter = calculate_ligand_diameter(args.ligand)
            # Set box size to sqrt(2) times the diameter
            box_dim = ligand_diameter * math.sqrt(2)
            # Add a small buffer of 2 Angstroms to be safe
            box_dim += 2.0
            print('[ftlk] LigandB max diameter calculated: %.3f Angstroms' % ligand_diameter)
            print('[ftlk] Auto-calculated box size (diameter * sqrt(2) + 2 A buffer): %.3f Angstroms' % box_dim)
            
        box = (box_dim, box_dim, box_dim)
        print('[ftlk] Box centered on original anchor atom %d.' % new_atomsA_serials[0])
        print('[ftlk] Box center: %.3f %.3f %.3f' % center)
        print('[ftlk] Box size  : %.3f %.3f %.3f' % box)

        print('[ftlk] Starting AGFR workflow to build target file...')
        trg_cov_base = os.path.join(tmpdir, 'target_cov')
        chain_id, res_name, res_seq = find_residue_info(final_rec_file, new_atomsA_serials[0])
        cov_residue_str = '%s:%s%s' % (chain_id, res_name, res_seq)
        
        agfr_cov_cmd = [
            PYTHON_EXE, os.path.join(PKG_BIN, 'runAGFR.py'), 
            '-r', final_rec_file, 
            '-b', 'user', '%.3f' % center[0], '%.3f' % center[1], '%.3f' % center[2],
            '%.3f' % box[0], '%.3f' % box[1], '%.3f' % box[2], 
            '-o', trg_cov_base,
            '-c', str(new_atomsA_serials[0]), str(new_atomsA_serials[1]),
            '-t', str(new_atomsA_serials[2]),
            '-x', cov_residue_str
        ]
        subprocess.check_call(agfr_cov_cmd)
        trg_cov_file = trg_cov_base + '.trg'
        print('[ftlk] Covalent target created: %s' % trg_cov_file)

        # --- FIX FOR jobName BEHAVIOR ---
        # Default ligand path and search name
        adfr_ligand_path = ligandB_with_sidechain_path
        base_name_for_search = os.path.basename(adfr_ligand_path).replace('.pdbqt', '')

        if args.jobName:
            # If jobName is given, rename the temporary ligand to match it.
            # This forces runADFR.py's default naming scheme to produce the desired output filenames.
            renamed_ligand_path = os.path.join(tmpdir, args.jobName + '.pdbqt')
            os.rename(adfr_ligand_path, renamed_ligand_path)
            
            # Update variables to use the new name
            adfr_ligand_path = renamed_ligand_path
            base_name_for_search = args.jobName
            print('[ftlk] Renamed temporary ligand to %s to control output name.' % os.path.basename(adfr_ligand_path))
        # --- END OF FIX ---

        print('[ftlk] Starting final ADFR docking run...')
        adfr_cmd_list = [
            PYTHON_EXE, os.path.join(PKG_BIN, 'runADFR.py'),
            '-l', adfr_ligand_path, # Use the (potentially renamed) ligand path
            '-t', trg_cov_file,
            '-C'] + [str(x) for x in new_covalent_atom_serials_for_adfr] + \
           ['--nbRuns', str(args.nbRuns),
            '--maxEvals', str(args.maxEvals),
            '--seed', str(args.seed),
            '-O']
        
        # --- REMOVED: Do NOT pass --jobName to runADFR.py anymore ---
        # if args.jobName:
        #     adfr_cmd_list.extend(['--jobName', args.jobName])
        
        adfr_cmd = adfr_cmd_list + extra
        
        print('[ftlk] >>', ' '.join(adfr_cmd))
        subprocess.check_call(adfr_cmd, cwd=output_dir)

        print('\n[ftlk] Docking finished.')

        # Use the base_name_for_search variable determined before the docking call
        search_dir = output_dir if output_dir else '.'
        search_pattern = os.path.join(search_dir, base_name_for_search + "*_out.pdbqt")
        
        print('[ftlk] Searching for result file with pattern: %s' % search_pattern)
        found_files = glob.glob(search_pattern)
        
        docking_output_file = None
        if found_files:
            docking_output_file = found_files[0]
            print('[ftlk] Found docking result: %s' % docking_output_file)

            # --- FIX for final filenames ---
            # Use the clean base name for final outputs
            final_cleaned_output = os.path.join(search_dir, "CLEANED_" + base_name_for_search + ".pdbqt")
            run_postprocessing(docking_output_file, final_cleaned_output)

            final_complex_output_pdb = os.path.join(search_dir, "COMPLEX_" + base_name_for_search + ".pdb")
            # --- END OF FIX ---

            write_combined_complex(
                args.complex, 
                docking_output_file, 
                final_complex_output_pdb,
                args.atomsA,
                chain_id,
                res_name,
                res_seq
            )
        else:
            print('[ftlk] WARNING: Docking output file matching pattern "%s" not found. Skipping final cleanup step.' % search_pattern)

    except Exception:
        print('\n[ftlk] ERROR: Docking aborted due to exception:')
        traceback.print_exc()
        if not args.keepTmp:
            args.keepTmp = True
            print('[ftlk] Temporary directory will be preserved for debugging.')
        print('[ftlk] Temporary directory has been saved to:', tmpdir)
        sys.exit(1)

    finally:
        if not args.keepTmp and 'tmpdir' in locals() and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == '__main__':
    main()