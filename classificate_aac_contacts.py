# Check possible contacts between structure and ligand
# python script.py ref.pdb LIG PDB_LIG.csv CHAIN(optional)
# GitHub: github.com/sulfierry/

from Bio.PDB import *
import numpy as np
import csv
import sys


# Variaveis globais ################################################################################################

# Distance threshold for hydrophobic interactions
hydrophobic_distance_threshold = 4.0

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
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
    
    # Átomos de carbono em grupos alquila e anéis alifáticos
    "C",   
    # Nome comum para carbonos em anéis aromáticos
    "CA", 
    # Outros nomes para carbonos em diferentes ambientes 
    "CH",
    # Carbonos numerados, comuns em moléculas pequenas 
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
    # Flúor
    "F",
    # Cloro   
    "Cl",
    # Bromo 
    "Br",
    # Iodo  
    "I"    
]


# Common Lennard-Jones parameters and charges for atom types from MMFF94 force field
atom_parameters = {
    'C': {'sigma': 3.50, 'epsilon': 0.066, 'charge': 0.0},   # Aproximado para carbono sp3
    'O': {'sigma': 3.07, 'epsilon': 0.152, 'charge': -0.5}, # Aproximado para oxigênio sp3
    'N': {'sigma': 3.25, 'epsilon': 0.170, 'charge': -0.5}, # Aproximado para nitrogênio sp3
    'H': {'sigma': 2.42, 'epsilon': 0.03, 'charge': 0.3},   # Hidrogênio
    'S': {'sigma': 3.80, 'epsilon': 0.250, 'charge': 0.0},  # Enxofre sp3
    'P': {'sigma': 3.74, 'epsilon': 0.200, 'charge': 0.5},  # Fósforo
    'F': {'sigma': 2.94, 'epsilon': 0.061, 'charge': -0.8}, # Flúor
    'Cl': {'sigma': 3.40, 'epsilon': 0.276, 'charge': -0.7},# Cloro
    'Br': {'sigma': 3.80, 'epsilon': 0.389, 'charge': -0.7},# Bromo
    'I':  {'sigma': 4.17, 'epsilon': 0.468, 'charge': -0.4} # Iodo
    # Adicione mais tipos atômicos conforme necessário
}

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

hydrogen_bond_acceptors = ["O", "N", "F"]
    
ionic_interactions = {

        "NA": ["F", "Cl", "Br", "I", "O"],
        "MG": ["F", "Cl", "Br", "I", "O"],
        "K" : ["F", "Cl", "Br", "I", "O"],
        "CA": ["F", "Cl", "Br", "I", "O"],
}

 ######################################################################################################################


def find_molecule(input_pdb, input_molecule):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("name", input_pdb)

    # Check if input_molecule is in molecule_class. If not, add it.
    if input_molecule not in molecule_class:
        molecule_class[input_molecule] = input_molecule

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
                if residue.get_resname() == input_molecule:
                    ligand_residue = residue
                    chain_select = chain.get_id()
                    break
            if ligand_residue:
                break

    if ligand_residue is None:
        raise ValueError(f"Ligand {input_molecule} not found in structure")
    
    return ligand_residue



def lennard_jones_potential(epsilon, sigma, r):
    """Calculate Lennard-Jones potential between two atoms."""
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)



# Define a function to check for potential hydrogen bonds
def is_interaction(atom1, atom2, residue_name, distance):

    if distance < 3.0:
        # Check for ionic interaction
        if residue_name in ionic_interactions:
            if atom1.get_name().startswith(tuple(
                ionic_interactions[residue_name])) or atom2.get_name().startswith(tuple(ionic_interactions[residue_name])):
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

    # Calculate Lennard-Jones potentials
    v_lj = lennard_jones_potential(
        # Este cálculo refere-se à combinação dos parâmetros de profundidade 
        # do poço de energia epsilon para dois átomos diferentes quando se modela 
        # uma interação via potencial de Lennard-Jones. A combinação geométrica (média geométrica) 
        # é comum para este parâmetro.
        (params1.get('epsilon', 0) * params2.get('epsilon', 0))**0.5,

        # Este cálculo refere-se à combinação dos parâmetros de distância de sigma para os mesmos dois átomos. 
        # sigma é geralmente interpretado como a distância em que o potencial interatômico entre dois átomos neutros é zero.
        # A combinação aritmética (média aritmética) é típica para este parâmetro.
        (params1.get('sigma', 0) + params2.get('sigma', 0)) / 2,
        distance
    )

    # Check for potential van der Waals interaction based on Lennard-Jones potential
    
    """
    Valor Conservador: Se você deseja ser mais conservador e focar apenas nas interações 
    mais fortes de van der Waals, pode considerar um valor limiar de −0.5 kcal/mol ou mais negativo.
    
    Valor Moderado: Um valor de −0.2 a −0.3 kcal/mol pode ser uma abordagem intermediária, 
    onde você identifica interações que têm uma contribuição notável, mas não são extremamente fracas.

    Análise Detalhada: Se o objetivo é uma análise mais detalhada e abrangente das interações, 
    incluindo as mais fracas, então −0.1 kcal/mol ou até um pouco mais positivo pode ser aceitável. 
    No entanto, essas interações devem ser interpretadas com cautela e corroboradas com outras evidências ou análises.
    
    """
    if v_lj < -0.1:  # using -0.1 kcal/mol as a threshold
        return "van der Waals"
    
    return "Non-specific"


def verify_near_residues(input_pdb, ligand_residue, treshold_distance):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("name", input_pdb)

    # List to store the next residues and the interaction
    near_residues = []

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

    return near_residues

def set_output(output_name, near_residues):

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
    
    return ("Successfully processed and saved!")



if __name__ == "__main__":


    input_pdb = sys.argv[1]
    input_molecule = str(sys.argv[2])
    output_name = str(sys.argv[3])

    # Distance from the select molecule
    treshold_distance = 4.0 

    ligand_residue = find_molecule(input_pdb, input_molecule)
    near_residues  = verify_near_residues(input_pdb, ligand_residue, treshold_distance)
    save_csv = set_output(output_name, near_residues)
