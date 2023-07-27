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
    "HOH": "polar tendo H carga parcial positiva, e O carga parcial negativa",
    "ACP": "ATP",
    "TPS": "TMP"
}

# Selecione o ligante
ligand_residue = None
for chain in structure[0]:
    for residue in chain:
        if residue.get_resname() == 'TPS':
            ligand_residue = residue
            break
    if ligand_residue:
        break

# Lista para armazenar os resíduos próximos
near_residues = []

# Verifica todos os átomos de todos os resíduos
for chain in structure[0]:
    for residue in chain:
        if residue != ligand_residue:  # Ignora o ligante
            for atom in residue:
                for ligand_atom in ligand_residue:
                    distance = atom - ligand_atom
                    if distance <= 4.0:
                        if residue not in near_residues:
                            near_residues.append(residue)
                        break  # Não precisa verificar os outros átomos deste resíduo

# Escreve os resíduos próximos em um arquivo csv
with open('output.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Nome do aminoácido", "Número", "Classificação"])
    for residue in near_residues:
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        aa_class = molecule_class.get(aa_name, "desconhecido")
        writer.writerow([aa_name, aa_num, aa_class])
        print("{:<20} {:<10} {:<30}".format(aa_name, aa_num, aa_class))
