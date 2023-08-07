# Check possible contacts between structure and ligand
# python script.py ref.pdb LIG PDB_LIG.csv CHAIN(optional)
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
    "K"  : "cofactor charge: +",
    "HOH": "polar charge: +-",
    "ACP": "ATP",
    "TPS": "TMP",
    "LIG": "LIG",
    "ANP": "ATP"
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
    raise ValueError(f"Ligand {molecule_select} not found in structure")


# List to store the next residues and the interaction
near_residues = []
# List to store the positions of the next residues
near_residues_positions = set()  


# Check all atoms of all residues
for chain in structure[0]:
    for residue in chain:
        if residue != ligand_residue:  # Ignore ligand
            min_distance = np.inf
            interacting_atoms = None
            for atom in residue:
                for ligand_atom in ligand_residue:
                    distance = atom - ligand_atom
                    if distance < min_distance:
                        min_distance = distance
                        interacting_atoms = (atom, ligand_atom)

            if min_distance <= 4.0:
                position = residue.get_full_id()
                if position not in [r[0].get_full_id() for r in near_residues]:
                    near_residues.append((residue, min_distance, interacting_atoms))

# Sort the near_residues list based on the distance (second element of the tuple)
near_residues.sort(key=lambda x: x[1])

# Write the residuals to a csv file
with open(output_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Amino acid", "Number", "Chain", "Classification", "Nearby atoms", "Distance Å"])  # Add "Chain"
    columns = ["Amino acid", "Number", "Chain", "Classification", "Nearby atoms", "Distance Å"]  # Add "Chain"
    print("{:^20} {:^10} {:^5} {:^30} {:^20} {:^10}".format(*columns))   # Add another column for "Chain"
    for residue, distance, atoms in near_residues:
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        chain_id = residue.get_parent().get_id()  # Get chain ID
        aa_class = molecule_class.get(aa_name, "unknown")
        atom1 = atoms[0].get_name() + "(" + residue.get_resname() + ")"
        atom2 = atoms[1].get_name() + "(" + ligand_residue.get_resname() + ")"
        nearby_atoms_str = "{:<10}-{:>10}".format(atom1, atom2)
        writer.writerow([aa_name, aa_num, chain_id, aa_class, nearby_atoms_str, round(distance, 2)])  # Add chain_id
        print("{:^20} {:^10} {:^5} {:^30} {:^20} {:^10.2f}".format(aa_name, aa_num, chain_id, aa_class, nearby_atoms_str, distance))  # Add chain_id
