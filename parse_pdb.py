def parse_pdb_extended_v5(pdb_file):
    """
    Parses a PDB file to extract chains, cofactors, ligands, and atom details.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: A dictionary containing lists:
              - 'chains': List of ATOM details.
              - 'cofactors': List of HETATM details matching cofactor names.
              - 'ligands': List of HETATM details matching ligand names.
    """
    
    def extract_atom_details(line):
        """Helper function to extract atom details from a line."""
        return {
            'serial_number': int(line[6:11].strip()),
            'name': line[12:16].strip(),
            'alt_loc': line[16].strip(),
            'res_name': line[17:20].strip(),
            'chain_id': line[21].strip(),
            'res_seq': int(line[22:26].strip()),
            'icode': line[26].strip(),
            'coord': [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ],
            'occupancy': float(line[54:60].strip()),
            'temp_factor': float(line[60:66].strip()),
            'element': line[76:78].strip(),
            'charge': line[78:80].strip()
        }

    chains = []
    cofactors = []
    ligands = []

    cofactor_names = ["MG", "ZN", "CA", "K", "NA"]
    ligand_names = ["ACP", "TPS"]

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chains.append(extract_atom_details(line))
            elif line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                details = extract_atom_details(line)
                if residue_name in cofactor_names:
                    cofactors.append(details)
                elif residue_name in ligand_names:
                    ligands.append(details)

    return {
        'chains': chains,
        'cofactors': cofactors,
        'ligands': ligands
    }


def print_pdb_structure(pdb_dict):
    """
    Prints the PDB structure from the dictionary in an organized manner.

    Parameters:
        pdb_dict (dict): Dictionary containing the parsed PDB lists.

    Returns:
        None: Simply prints the PDB structured data.
    """

def format_line(atom_data, atom_type="ATOM  "):
    return (
        f"{atom_type:6s}{atom_data['serial_number']:5d} {atom_data['name']:<4s} {atom_data['alt_loc']:1s}{atom_data['res_name']:<3s} "
        f"{atom_data['chain_id']:1s}{atom_data['res_seq']:4d}{atom_data['icode']:1s}   "
        f"{atom_data['coord'][0]:8.3f}{atom_data['coord'][1]:8.3f}{atom_data['coord'][2]:8.3f}"
        f"{atom_data['occupancy']:6.2f}{atom_data['temp_factor']:6.2f}          "
        f"{atom_data['element']:^2s}{atom_data['charge']:2s}\n"
    )

def print_pdb_structure(pdb_dict):
    # Print the chains first
    for atom in pdb_dict['chains']:
        print(format_line(atom), end='')

    # Print cofactors
    print("TER")
    for atom in pdb_dict['cofactors']:
        print(format_line(atom, "HETATM"), end='')

    # Print ligands
    print("TER")
    for atom in pdb_dict['ligands']:
        print(format_line(atom, "HETATM"), end='')

    print("END")



# Parse the provided PDB
pdb_dict = parse_pdb_extended_v5("../model_structure/3c9t.pdb")
print_pdb_structure(pdb_dict)
#print(type(pdb_dict))
