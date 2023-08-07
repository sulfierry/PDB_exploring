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


# Amino acids and nucleic acid bases classification
molecule_class = {
    # Amino Acids
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
    
    # Nucleic Acid Bases
    "A": "adenine base",
    "C": "cytosine base",
    "G": "guanine base",
    "T": "thymine base",   # DNA
    "U": "uracil base",    # RNA

    # Cofactors
    "Mg" : "cofactor charge: +",
    "K"  : "cofactor charge: +",
    "Ca" : "cofactor charge +",
    "Na" : "cofactor charge",
    "HoH": "partial polar charge"
}


# Check if molecule_to_check is in molecule_class. If not, add it.
if molecule_select not in molecule_class:
    molecule_class[molecule_select] = molecule_select

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

# Define a function to check for potential hydrogen bonds
def is_interaction(atom1, atom2, residue_name, distance):
    hydrogen_bond_acceptors = ["O", "N", "F"]
    ionic_interactions = {
        "NA": ["F", "Cl", "Br", "I", "O", "S"],
        "MG": ["F", "Cl", "Br", "I", "O", "S"],
        "K": ["F", "Cl", "Br", "I", "O", "S"],
        "CA": ["F", "Cl", "Br", "I", "O", "S"],
        #"HOH": ["F", "Cl", "Br", "I", "O", "S"]
    }

    if distance < 3.0:
        # Check for ionic interaction
        if residue_name in ionic_interactions:
            if atom1.get_name().startswith(tuple(ionic_interactions[residue_name])) or atom2.get_name().startswith(tuple(ionic_interactions[residue_name])):
                return "Ionic"

        # Check for hydrogen bond
        if atom1.get_name().startswith(tuple(hydrogen_bond_acceptors)) and atom2.get_name().startswith(tuple(hydrogen_bond_acceptors)):
            return "True"
    
    return "False"

# ...




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
    # Change the order of columns here
    writer.writerow(["Molecule", "Classification", "Number", 
                     "Chain", "Nearby atoms", "Distance(Å)", "H bond"]) 
    
    columns = ["Molecule", "Classification", "Number", 
               "Chain", "Nearby atoms", "Distance(Å)", "H bond"] 
    
    print("{:^20} {:^30} {:^10} {:^5} {:^20} {:^10} {:^6}".format(*columns))
    for residue, distance, atoms in near_residues:
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        chain_id = residue.get_parent().get_id()  # Get chain ID

        if aa_name in ["NA", "MG", "K", "CA", "HOH"]:
            aa_class = "cofactor"
        else:
            aa_class = molecule_class.get(aa_name, "unknown")

        atom1 = atoms[0].get_name() + "(" + residue.get_resname() + ")"
        atom2 = atoms[1].get_name() + "(" + ligand_residue.get_resname() + ")"
        nearby_atoms_str = "{:<10}-{:>10}".format(atom1, atom2)
        
        # Check for hydrogen bond
        h_bond = is_interaction(atoms[0], atoms[1], residue.get_resname(), distance)

        # If the result is "probable", set h_bond to "probable". Otherwise, set it to "True" or "False".
        #h_bond = "probable" if h_bond_check == "probable" else ("True" if h_bond_check else "False")

        # Change the order of written rows here to match column header order
        writer.writerow([aa_name, aa_class, aa_num, 
                         chain_id, nearby_atoms_str, round(distance, 2), h_bond]) 
        
        print("{:^20} {:^30} {:^10} {:^5} {:^20} {:^10.2f} {:^6}".format(
            aa_name, aa_class, aa_num, chain_id, nearby_atoms_str, distance, h_bond))
