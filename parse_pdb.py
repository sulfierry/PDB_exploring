
def distance(atom1, atom2):
    return ((atom1['x'] - atom2['x'])**2 + (atom1['y'] - atom2['y'])**2 + (atom1['z'] - atom2['z'])**2)**0.5

# Parser for PDB structure with consideration for protein structure, ligands, and cofactors
def parse_pdb_extended_v2(pdb_path):
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    structure = {'chains': {}, 'cofactors': [], 'ligand': {'resname': 'LIG', 'atoms': []}}
    cofactor_names = ["MG", "NA", "CA", "K", "FE"]
    waters_near_cofactors = []

    for line in lines:
        if line.startswith("ATOM"):
            chain_id = line[21]
            residue_num = int(line[22:26].strip())
            residue_name = line[17:20].strip()
            if chain_id not in structure['chains']:
                structure['chains'][chain_id] = {'residues': {}}
            if residue_num not in structure['chains'][chain_id]['residues']:
                structure['chains'][chain_id]['residues'][residue_num] = {'resname': residue_name, 'atoms': []}
            atom_data = {
                'atom_name': line[12:16].strip(),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54])
            }
            structure['chains'][chain_id]['residues'][residue_num]['atoms'].append(atom_data)
        elif line.startswith("HETATM"):
            residue_name = line[17:20].strip()
            atom_data = {
                'atom_name': line[12:16].strip(),
                'chain_id': line[21],
                'residue_num': int(line[22:26].strip()),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54])
            }
            if residue_name in cofactor_names:
                cofactor_exists = False
                for cof in structure['cofactors']:
                    if cof['resname'] == residue_name and cof['chain_id'] == atom_data['chain_id'] and cof['residue_num'] == atom_data['residue_num']:
                        cof['atoms'].append(atom_data)
                        cofactor_exists = True
                        break
                if not cofactor_exists:
                    structure['cofactors'].append({
                        'resname': residue_name,
                        'chain_id': atom_data['chain_id'],
                        'residue_num': atom_data['residue_num'],
                        'atoms': [atom_data]
                    })
            elif residue_name == "HOH":
                waters_near_cofactors.append(atom_data)
            else:
                structure['ligand']['atoms'].append(atom_data)

    for water in waters_near_cofactors:
        for cofactor in structure['cofactors']:
            for atom in cofactor['atoms']:
                if distance(water, atom) <= 5.0:
                    cofactor_exists = False
                    for cof in structure['cofactors']:
                        if cof['resname'] == "HOH" and cof['chain_id'] == water['chain_id'] and cof['residue_num'] == water['residue_num']:
                            cof['atoms'].append(water)
                            cofactor_exists = True
                            break
                    if not cofactor_exists:
                        structure['cofactors'].append({
                            'resname': "HOH",
                            'chain_id': water['chain_id'],
                            'residue_num': water['residue_num'],
                            'atoms': [water]
                        })
                    break
    return structure



def print_chains(structure_dict):
    atom_count = 1
    for chain_id, chain_data in structure_dict['chains'].items():
        for residue_num, residue_data in chain_data['residues'].items():
            for atom in residue_data['atoms']:
                print(f"ATOM  {atom_count:5} {atom['atom_name']:^4} {residue_data['resname']:<3} {chain_id} {residue_num:4}    {atom['x']:8.3f} {atom['y']:8.3f} {atom['z']:8.3f}")
                atom_count += 1

def print_cofactors(structure_dict):
    atom_count = 1
    for cofactor in structure_dict['cofactors']:
        for atom in cofactor['atoms']:
            print(f"HETATM{atom_count:5} {atom['atom_name']:^4} {cofactor['resname']:<3} {atom['chain_id']} {atom['residue_num']:4}    {atom['x']:8.3f} {atom['y']:8.3f} {atom['z']:8.3f}")
            atom_count += 1

def print_ligand(structure_dict):
    atom_count = 1
    ligand = structure_dict['ligand']
    for atom in ligand['atoms']:
        print(f"HETATM{atom_count:5} {atom['atom_name']:^4} {ligand['resname']:<3} {atom['chain_id']} {atom['residue_num']:4}    {atom['x']:8.3f} {atom['y']:8.3f} {atom['z']:8.3f}")
        atom_count += 1




# Parse the provided PDB
parsed_structure_extended = parse_pdb_extended_v2("../model_structure/3c9t.pdb")
print('Chains: \n') 
# Display the keys of the parsed structure for a preview
print_chains(parsed_structure_extended)
print('\n Cofactors: \n')
  # Display the keys of the parsed structure for a preview
print_cofactors(parsed_structure_extended)
print('\n Ligand: \n')
print_ligand(parsed_structure_extended)