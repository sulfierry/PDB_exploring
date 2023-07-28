# Check contacts between structure and ligand
# python script.py ref.pdb A LIG PDB_LIG.csv
# GitHub: github.com/sulfierry/

from Bio.PDB import *
import numpy as np
import csv
import sys

input_pdb = sys.argv[1]
molecule_select = str(sys.argv[2])
output_name = str(sys.argv[3])



parser = PDBParser(QUIET=True)
structure = parser.get_structure("name", input_pdb)


# Amino acids classification
molecule_class = {
    "ALA": "hydrophobic",
    "ILE": "hydrophobic",
    "LEU": "hydrophobic",
    "VAL": "hydrophobic",
    "PHE": "hydrophobic",
    "PRO": "hydrophobic",
    "TRP": "hydrophobic",
    "MET": "hydrophobic",
    "GLY": "hydrophobic",
    "CYS": "polar charge: 0",
    "SER": "polar charge: 0",
    "THR": "polar charge: 0",
    "TYR": "polar charge: 0",
    "ASN": "polar charge: 0",
    "GLN": "polar charge: 0",
    "HIS": "polar charge: +",
    "LYS": "polar charge: +",
    "ARG": "polar charge: +",
    "ASP": "polar charge: -",
    "GLU": "polar charge: -",
    "MG" : "cofactor charge: +",
    "HOH": "polar charge: +-",
    "ACP": "ATP",
    "TPS": "TMP"
}

# Ligand selection
# a tuple with ('string', residue number)
ligand_residue = None
chain_select = None

# Check if the chain_select was provided as an argument
if len(sys.argv) > 4:
    chain_select = str(sys.argv[4])

# Look for the ligand in the specified chain or in all chains
for model in structure:
    for chain in model:
        # If chain_select was provided, check only that chain
        if chain_select and chain_select.upper() != chain.get_id():
            continue

        for residue in chain:
            if residue.get_resname() == molecule_select:
                ligand_residue = residue
                chain_select = chain.get_id()
                break
        if ligand_residue:
            break

if ligand_residue is None:
    raise ValueError(f"Ligante {molecule_select} não encontrado na estrutura")


# List to store the next residues and the interaction
near_residues = []
# List to store the positions of the next residues
near_residues_positions = set()  


# Check all atoms of all residues
for chain in structure[0]:
    for residue in chain:
        if residue != ligand_residue:  
            # Ignore ligand
            min_distance = np.inf
            interacting_atoms = None
            for atom in residue:
                for ligand_atom in ligand_residue:
                    distance = atom - ligand_atom
                    if distance < min_distance:
                        # Keep the minimum distance
                        min_distance = distance  
                        # Keeps the atoms that produced the minimum distance
                        interacting_atoms = (atom, ligand_atom)  

            if min_distance <= 4.0:
                interaction = 'hydrogen bonding: ' + str(round(min_distance,2)) + ' Å' if min_distance <= 2.4 else 'weak hydrogen bonding: ' + str(round(min_distance,2)) + ' Å' if min_distance <= 2.9 else 'van der Waals: ' + str(round(min_distance,2)) + ' Å'
                hydrophobic_interaction = "Yep (" + str(round(min_distance,2)) + ' Å)' if min_distance <= 4.0 and molecule_class[ligand_residue.get_resname()] == "hydrophobic" else "Not "
                # Gets the unique position of the residue in the structure
                position = residue.get_full_id()  
                if position not in [r[0].get_full_id() for r in near_residues]:
                    near_residues.append((residue, interaction, interacting_atoms, hydrophobic_interaction))


# Write the next residuals to a csv file
with open(output_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Amino acid", "Number", "Classification", 
                     "Intermolecular interaction", "hydrophobic interaction", "Interacting atoms"])
    columns = ["Amino acid", "Number", "Classification", 
               "Intermolecular interaction", "hydrophobic interaction", "Interacting atoms"]
    print("{:^20} {:^10} {:^30} {:^40} {:^20} {:^20}".format(*columns))   
    for residue, interaction, atoms, hydrophobic_interaction in near_residues:
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        aa_class = molecule_class.get(aa_name, "unknown")
        atom1 = atoms[0].get_name() + "(" + residue.get_resname() + ")"
        atom2 = atoms[1].get_name() + "(" + ligand_residue.get_resname() + ")"
        interacting_atoms_str = "{:<10}-{:>10}".format(atom1, atom2)  # Adjust the size of the atoms fields
        writer.writerow([aa_name, aa_num, aa_class, interaction, hydrophobic_interaction, interacting_atoms_str])
        print("{:^20} {:^10} {:^30} {:^40} {:^22} {:^25}".format(aa_name, aa_num, aa_class, interaction, hydrophobic_interaction, interacting_atoms_str))

