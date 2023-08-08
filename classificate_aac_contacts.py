# Check possible contacts between structure and ligand
# python script.py ref.pdb LIG PDB_LIG.csv CHAIN(optional)
# GitHub: github.com/sulfierry/

from Bio.PDB import *
import numpy as np
import csv
import sys


# Hydrophobic residues and atoms for identifying hydrophobic interactions
hydrophobic_residues = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "TYR"]
hydrophobic_atoms = [
    # Alanina
    "CB",
    # Valina
    "CB", "CG1", "CG2",
    # Isoleucina
    "CB", "CG1", "CG2", "CD1",
    # Leucina
    "CB", "CG", "CD1", "CD2",
    # Metionina
    "CB", "CG", "SD", "CE",
    # Prolina
    "CB", "CG", "CD",
    # Fenilalanina
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
    # Triptofano
    "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2",
    # Tirosina
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH",

    "C",   # Átomos de carbono em grupos alquila e anéis alifáticos
    "CA",  # Nome comum para carbonos em anéis aromáticos
    "CH",  # Outros nomes para carbonos em diferentes ambientes
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",  # Carbonos numerados, comuns em moléculas pequenas
    "F",   # Flúor
    "Cl",  # Cloro
    "Br",  # Bromo
    "I"    # Iodo
]

# Distance from the select molecule
treshold_distance = 4.0 
# Distance threshold for hydrophobic interactions
hydrophobic_distance_threshold = 4.0


input_pdb = sys.argv[1]
molecule_select = str(sys.argv[2])
output_name = str(sys.argv[3])


parser = PDBParser(QUIET=True)
structure = parser.get_structure("name", input_pdb)



# Constants
epsilon_0 = 8.854187817e-12  # permittivity of vacuum in C^2/N*m^2
four_pi_epsilon_0 = 4 * 3.14159265359 * epsilon_0  # 4*pi*epsilon_0 for Coulomb's law

# Common Lennard-Jones parameters and charges for atom types from CHARMM force field
atom_parameters = {
    'C': {'sigma': 3.55, 'epsilon': 0.066, 'charge': 0.0},
    'O': {'sigma': 3.12, 'epsilon': 0.12, 'charge': -0.5},
    'N': {'sigma': 3.25, 'epsilon': 0.17, 'charge': -0.5},
    'H': {'sigma': 2.42, 'epsilon': 0.03, 'charge': 0.3},
    'S': {'sigma': 3.79, 'epsilon': 0.25, 'charge': 0.0},
    'P': {'sigma': 3.74, 'epsilon': 0.20, 'charge': 0.5},
    # Add more atom types as necessary
}

def coulomb_potential(q1, q2, r):
    """Calculate Coulomb potential between two charges."""
    return (q1 * q2) / (four_pi_epsilon_0 * r)

def lennard_jones_potential(epsilon, sigma, r):
    """Calculate Lennard-Jones potential between two atoms."""
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

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
    "HOH": "partial polar charge"
}


# Check if molecule_select is in molecule_class. If not, add it.
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
        "NA": ["F", "Cl", "Br", "I", "O"],
        "MG": ["F", "Cl", "Br", "I", "O"],
        "K" : ["F", "Cl", "Br", "I", "O"],
        "CA": ["F", "Cl", "Br", "I", "O"],
    }

    if distance < 3.0:
        # Check for ionic interaction
        if residue_name in ionic_interactions:
            if atom1.get_name().startswith(tuple(ionic_interactions[residue_name])) or atom2.get_name().startswith(tuple(ionic_interactions[residue_name])):
                return "Ionic"

        # Check for hydrogen bond
        if atom1.get_name().startswith(tuple(hydrogen_bond_acceptors)) and atom2.get_name().startswith(tuple(hydrogen_bond_acceptors)):
            return "Hydrogen bond"
            
    
    # Check for hydrophobic interactions
    if residue_name in hydrophobic_residues:
        if atom1.get_name() in hydrophobic_atoms or atom2.get_name() in hydrophobic_atoms:
            if distance <= hydrophobic_distance_threshold:
                return "Hydrophobic"

    # Check for van der Waals interaction
    atom1_type = atom1.get_name()[0]  # Simplified to get first letter, might need refinement
    atom2_type = atom2.get_name()[0]

    # Get atom parameters
    params1 = atom_parameters.get(atom1_type, {})
    params2 = atom_parameters.get(atom2_type, {})

    # Calculate Coulomb and Lennard-Jones potentials
    v_coulomb = coulomb_potential(params1.get('charge', 0), params2.get('charge', 0), distance)
    v_lj = lennard_jones_potential(
        (params1.get('epsilon', 0) * params2.get('epsilon', 0))**0.5,
        (params1.get('sigma', 0) + params2.get('sigma', 0)) / 2,
        distance
    )

    # Check for potential van der Waals interaction based on Lennard-Jones potential
    if v_lj < -0.1:  # using -0.1 kcal/mol as a threshold
        return "van der Waals"
    
    return "Non-specific"


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

            if min_distance <= treshold_distance:
                position = residue.get_full_id()
                if position not in [r[0].get_full_id() for r in near_residues]:
                    near_residues.append((residue, min_distance, interacting_atoms))

# Sort the near_residues list based on the distance (second element of the tuple)
near_residues.sort(key=lambda x: x[1])


with open(output_name, 'w', newline='') as file:

    writer = csv.writer(file)
    # Change the order of columns here
    writer.writerow(["Molecule", "Classification", "Number", 
                     "Chain", "Nearby atoms", "Distance(Å)", "Interaction"]) 
    
    columns = ["Molecule", "Classification", "Number", 
               "Chain", "Nearby atoms", "Distance(Å)", "Interaction"] 
    
    print("{:^20} {:^30} {:^10} {:^5} {:^20} {:^10} {:^20}".format(*columns))

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
        probable_interaction = is_interaction(atoms[0], atoms[1], residue.get_resname(), distance)

        # Change the order of written rows here to match column header order
        writer.writerow([aa_name, aa_class, aa_num, 
                         chain_id, nearby_atoms_str, round(distance, 2), probable_interaction]) 
        
        print("{:^20} {:^30} {:^10} {:^5} {:^20} {:^10.2f} {:^20}".format(
            aa_name, aa_class, aa_num, chain_id, nearby_atoms_str, distance, probable_interaction))
