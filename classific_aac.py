from Bio.PDB import *
import numpy as np
import csv

# Abre o arquivo PDB e pega a estrutura
parser = PDBParser(QUIET=True)
structure = parser.get_structure("name", "./3c9t.pdb")

# Classificação dos aminoácidos
molecule_class = {
    "ALA": "hidrofóbico",
    "ILE": "hidrofóbico",
    "LEU": "hidrofóbico",
    "VAL": "hidrofóbico",
    "PHE": "hidrofóbico",
    "PRO": "hidrofóbico",
    "TRP": "hidrofóbico",
    "MET": "hidrofóbico",
    "GLY": "hidrofóbico",
    "CYS": "polar sem carga",
    "SER": "polar sem carga",
    "THR": "polar sem carga",
    "TYR": "polar sem carga",
    "ASN": "polar sem carga",
    "GLN": "polar sem carga",
    "HIS": "polar com carga positiva",
    "LYS": "polar com carga positiva",
    "ARG": "polar com carga positiva",
    "ASP": "polar com carga negativa",
    "GLU": "polar com carga negativa",
    "MG" : "cofator com carga positiva",
    "HOH": "polar carga parcial negativa",
    "ACP": "ATP",
    "TPS": "TMP"
}

# Selecione o ligante
ligand_selection = ('A', 'TPS')  # uma tupla com ('cadeia', número do resíduo)
ligand_residue = None

for chain in structure[0]:
    if chain.get_id() == ligand_selection[0]:  # checa a cadeia
        for residue in chain:
            if isinstance(ligand_selection[1], str) and residue.get_resname() == ligand_selection[1]:
                ligand_residue = residue
                break
            elif isinstance(ligand_selection[1], int) and residue.get_id()[1] == ligand_selection[1]:
                ligand_residue = residue
                break
        if ligand_residue:
            break

if ligand_residue is None:
    raise ValueError(f"Ligante {ligand_selection} não encontrado na estrutura")


# Dicionário para armazenar os resíduos próximos e a interação
near_residues = {}

# Verifica todos os átomos de todos os resíduos
for chain in structure[0]:
    for residue in chain:
        if residue != ligand_residue:  # Ignora o ligante
            min_distance = np.inf
            interacting_atoms = None
            for atom in residue:
                for ligand_atom in ligand_residue:
                    distance = atom - ligand_atom
                    if distance < min_distance:
                        min_distance = distance  # Mantém a distância mínima
                        interacting_atoms = (atom, ligand_atom)  # Mantém os átomos que produziram a distância mínima
            if min_distance <= 4.0:
                interaction = 'H-bond (1.5 to 2.7 Å)' if min_distance <= 2.9 else 'van der Waals (3.0 to 4.0 Å)'
                if residue not in near_residues or near_residues[residue][0] == 'van der Waals (3.0 to 4.0 Å)':
                    near_residues[residue] = (interaction, interacting_atoms)

# Escreve os resíduos próximos em um arquivo csv
with open('output.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Nome do aminoácido", "Número", "Classificação", "Interação", "Átomos Interagindo"])
    print("{:<20} {:<10} {:<30} {:<15} {:<20}".format("Aminoácido", "Número", "Classificação", "Interação", "Átomos Interagindo"))   
    for residue, (interaction, atoms) in near_residues.items():
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        aa_class = molecule_class.get(aa_name, "desconhecido")
        atom1 = atoms[0].get_name() + "(" + residue.get_resname() + ")"
        atom2 = atoms[1].get_name() + "(" + ligand_residue.get_resname() + ")"
        interacting_atoms_str = "{:<10}-{:>10}".format(atom1, atom2)  # Ajusta o tamanho dos campos dos átomos
        writer.writerow([aa_name, aa_num, aa_class, interaction, interacting_atoms_str])
        print("{:<20} {:<10} {:<30} {:<15} {:<20}".format(aa_name, aa_num, aa_class, interaction, interacting_atoms_str))
