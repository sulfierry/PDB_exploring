import sys
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np

# Definindo os parâmetros de Lennard-Jones
# Nota: Esses valores são apenas placeholders, por favor, substitua-os pelos valores reais
lj_params = {
    'C': (3.55, 0.276),
    'H': (2.42, 0.062),
    'O': (3.12, 0.210),
    'N': (3.25, 0.170),
    # Adicione mais átomos conforme necessário
}

# Definindo as cargas parciais dos átomos
# Nota: Esses valores são apenas placeholders, por favor, substitua-os pelos valores reais
charges = {
    'C': -0.18,
    'H': 0.20,
    'O': -0.65,
    'N': -0.57,
    # Adicione mais átomos conforme necessário
}

ke = 332.06371  # Constante de Coulomb em kcal*mol^-1*Å*e^-2

# Demais partes do código...

# Checando todos os átomos de todos os resíduos
for chain in structure[0]:
    for residue in chain:
        if residue != ligand_residue:  # Ignora o ligante
            min_distance = np.inf
            interacting_atoms = None
            for atom in residue:
                for ligand_atom in ligand_residue:
                    distance = atom - ligand_atom
                    if distance < min_distance:
                        # Mantém a distância mínima
                        min_distance = distance
                        # Mantém os átomos que produziram a distância mínima
                        interacting_atoms = (atom, ligand_atom)

            if min_distance <= 4.0:
                interaction = None
                hydrophobic_interaction = "Not "
                if atom.get_name() in ['H', 'F', 'O', 'N']:
                    if min_distance <= 2.5:
                        interaction = 'hydrogen bonding: ' + str(round(min_distance,2)) + ' Å'
                    elif min_distance <= 2.9:
                        interaction = 'weak hydrogen bonding: ' + str(round(min_distance,2)) + ' Å'
                if interaction is None:  # Se a interação não foi definida acima, é de Van der Waals
                    # Calcula o potencial de Lennard-Jones
                    sigma1, epsilon1 = lj_params[atom.get_name()]
                    sigma2, epsilon2 = lj_params[ligand_atom.get_name()]
                    sigma = (sigma1 + sigma2) / 2
                    epsilon = np.sqrt(epsilon1 * epsilon2)
                    # Note: a distância precisa ser convertida para o mesmo sistema de unidades do sigma
                    interaction_energy = 4 * epsilon * ((sigma / min_distance)**12 - (sigma / min_distance)**6)
                    # Calcula a interação eletrostática de Coulomb
                    coulomb_energy = ke * charges[atom.get_name()] * charges[ligand_atom.get_name()] / min_distance
                    if interaction_energy < 0:  # Energia negativa significa atração
                        interaction = 'Van der Waals: ' + str(round(min_distance,2)) + ' Å, energy: ' + str(round(interaction_energy,2)) + ' kcal/mol, Coulomb: ' + str(round(coulomb_energy,2)) + ' kcal/mol'
                if min_distance <= 4.0 and molecule_class[ligand_residue.get_resname()] == "hydrophobic":
                    hydrophobic_interaction = "Yep (" + str(round(min_distance,2)) + ' Å)'
                # Obtém a posição única do resíduo na estrutura
                position = residue.get_full_id()
                if position not in [r[0].get_full_id() for r in near_residues]:
                    near_residues.append((residue, interaction, interacting_atoms, hydrophobic_interaction))
