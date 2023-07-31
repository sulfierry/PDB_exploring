import sys
import numpy as np
from Bio.PDB import *
import csv

input_pdb = sys.argv[1]
molecule_select = str(sys.argv[2])
output_name = str(sys.argv[3])

parser = PDBParser(QUIET=True)
structure = parser.get_structure("name", input_pdb)

# Amino acids classification
molecule_class = { ... }  # Inserir a classificação dos aminoácidos

# Ligand selection
# a tuple with ('string', residue number)
ligand_residue = None
chain_select = None

# Check if the chain_select was provided as an argument
if len(sys.argv) > 4:
    chain_select = str(sys.argv[4])

near_residues = []

# Lennard-Jones parameters
lj_params = {
    'C': (3.55, 0.105),   # Example values, replace with actual values
    'H': (2.42, 0.015),   # Example values, replace with actual values
    'O': (3.12, 0.16),    # Example values, replace with actual values
    'N': (3.25, 0.17),    # Example values, replace with actual values
    # Add more atoms as needed
}

# Atomic charges
charges = {
    'C': -0.18,  # Example values, replace with actual values
    'H': 0.20,   # Example values, replace with actual values
    'O': -0.65,  # Example values, replace with actual values
    'N': -0.57,  # Example values, replace with actual values
    # Add more atoms as needed
}

ke = 332.06371  # Coulomb constant in kcal*mol^-1*Å*e^-2

# Checking all atoms from all residues...
for chain in structure[0]:
    for residue in chain:
        if residue != ligand_residue:  # Ignore the ligand
            min_distance = np.inf
            interacting_atoms = None
            for atom in residue:
                for ligand_atom in ligand_residue:
                    distance = np.linalg.norm(np.array(atom.coord) - np.array(ligand_atom.coord))
                    if distance < min_distance:
                        # Keep minimum distance
                        min_distance = distance
                        # Keep the atoms that produced the minimum distance
                        interacting_atoms = (atom, ligand_atom)

            if min_distance <= 4.0:
                interaction = None
                hydrophobic_interaction = "Not "
                if atom.get_name() in lj_params and ligand_atom.get_name() in lj_params:
                    # Calculate van der Waals interaction using Lennard-Jones potential
                    sigma1, epsilon1 = lj_params[atom.get_name()]
                    sigma2, epsilon2 = lj_params[ligand_atom.get_name()]
                    sigma = (sigma1 + sigma2) / 2
                    epsilon = np.sqrt(epsilon1 * epsilon2)
                    vdw_energy = 4 * epsilon * ((sigma / min_distance)**12 - (sigma / min_distance)**6)
                    # Calculate electrostatic interaction using Coulomb's law
                    coulomb_energy = ke * charges[atom.get_name()] * charges[ligand_atom.get_name()] / min_distance
                    interaction = 'VDW: ' + str(round(vdw_energy,2)) + ' kcal/mol, Coulomb: ' + str(round(coulomb_energy,2)) + ' kcal/mol'
                
                if min_distance <= 4.0 and molecule_class[ligand_residue.get_resname()] == "hydrophobic":
                    hydrophobic_interaction = "Yep (" + str(round(min_distance,2)) + ' Å)'
                
                # Get the residue's unique position in the structure
                position = residue.get_full_id()
                if position not in [r[0].get_full_id() for r in near_residues]:
                    near_residues.append((residue, interaction, interacting_atoms, hydrophobic_interaction))

# Write the next residuals to a csv file
with open(output_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Amino acid", "Number", "Classification", "Probable Intermolecular interaction", "Hydrophobic interaction", "Interacting atoms"])
    columns = ["Amino acid", "Number", "Classification", "Probable intermolecular interaction", "Hydrophobic interaction", "Interacting atoms"]
    print("{:^20} {:^10} {:^30} {:^40} {:^20} {:^20}".format(*columns))   
    for residue, interaction, atoms, hydrophobic_interaction in near_residues:
        aa_name = residue.get_resname()
        aa_num = residue.get_id()[1]
        aa_class = molecule_class.get(aa_name, "unknown")
        atom1 = atoms[0].get_name() + "(" + residue.get_resname() + ")"
        atom2 = atoms[1].get_name() + "(" + ligand_residue.get_resname() + ")"
        interacting_atoms_str = "{:<10}-{:>10}".format(atom1, atom2)  # Adjust the size of the atoms fields
        writer.writerow([aa_name, aa_num, aa_class, interaction, hydrophobic_interaction, interacting_atoms
