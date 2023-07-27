from Bio.PDB import *
import numpy as np
import csv
import sys

input_pdb = sys.argv[1]
chain_select = str(sys.argv[2])
molecule_select = str(sys.argv[3])
output_name = str(sys.argv[4])


parser = PDBParser(QUIET=True)
structure = parser.get_structure("name", input_pdb)


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
    "HOH": "polar carga parcial neutra",
    "ACP": "ATP",
    "TPS": "TMP"
}

# Selecione o ligante
ligand_selection = (chain_select, molecule_select)  # uma tupla com ('cadeia', número do resíduo)
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
                interaction = 'H-bond (1.5 to 2.9 Å)' if min_distance <= 2.9 else 'van der Waals (3.0 to 4.0 Å)'
                hydrophobic_interaction = "Yes" if min_distance <= 4.0 and molecule_class[residue.get_resname()] == "hidrofóbico" else "No"
                if residue not in near_residues or near_residues[residue][0] == 'van der Waals (3.0 to 4.0 Å)':
                    near_residues[residue] = (interaction, interacting_atoms, hydrophobic_interaction)

# Escreve os resíduos próximos em um arquivo csv
with open(output_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Aminoácido", "Numero", "Classificação", "Interação intermolecular", "Interacção hidrofóbica", "Átomos interagindo"])
    print("{:<20} {:<10} {:<30} {:<30} {:<20} {:<20}".format("Aminoácido", "Número", "Classificação", "Interação intermolecular", "Interacção hidrofóbica", "Átomos interagindo"))   
    for residue, (interaction, atoms, hydrophobic_interaction) in near_residues.items():
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        aa_class = molecule_class.get(aa_name, "desconhecido")
        atom1 = atoms[0].get_name() + "(" + residue.get_resname() + ")"
        atom2 = atoms[1].get_name() + "(" + ligand_residue.get_resname() + ")"
        interacting_atoms_str = "{:<10}-{:>10}".format(atom1, atom2)  # Ajusta o tamanho dos campos dos átomos
        writer.writerow([aa_name, aa_num, aa_class, interaction, hydrophobic_interaction, interacting_atoms_str])
        print("{:<20} {:<10} {:<30} {:<30} {:<22} {:<25}".format(aa_name, aa_num, aa_class, interaction, hydrophobic_interaction, interacting_atoms_str))
