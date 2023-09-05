import numpy as np
import matplotlib.pyplot as plt

# Função de plot 3D para visualizar os voxels
def plot_voxels(voxel_grid):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.voxels(voxel_grid, edgecolor="k")
    plt.show()


def plot_voxel(voxel):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = np.where(voxel == 1)
    ax.scatter(x, y, z, alpha=0.6, s=15, c='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

# Parser for PDB structure with consideration for protein structure, ligands, and cofactors
def parse_pdb(pdb_file):
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
        try:
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
        except ValueError:
            print(f"Problem with line: {line}")
            raise


    chains = []
    cofactors = []
    ligands = []

    cofactor_names = ["MG", "ZN", "CA", "K", "NA", "FE", "CL", "HOH"]
    ligand_names = ["ACP", "TPS", "TMP", "TPP", "LIG"]

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

# Função para converter coordenadas do PDB para índices de voxel
def coord_to_voxel(coord, origin, grid_size):
    voxel_coord = np.floor((coord - origin) / grid_size).astype(int)
    if np.any(voxel_coord < 0) or np.any(voxel_coord >= grid_dim):
        print(f"Invalid voxel coordinates: {voxel_coord} for atom coordinate: {coord}")
        return None
    return voxel_coord

# Função para converter o PDB (no formato de dicionário) para um voxel grid
def pdb_to_voxel_with_info(parsed_pdb):
    voxel_grid = np.zeros(grid_dim, dtype= int)
    info_grid = np.empty(grid_dim, dtype=object)

    # Iterar sobre todos os átomos e converter suas coordenadas para coordenadas de voxel
    for atom_section in ["chains", "cofactors", "ligands"]:
        for atom_details in parsed_pdb[atom_section]:
            voxel_coord = coord_to_voxel(np.array(atom_details["coord"]), origin, grid_size)
            if voxel_coord is not None:
                voxel_grid[tuple(voxel_coord)] = 1
                info_grid[tuple(voxel_coord)] = atom_details

    return voxel_grid, info_grid



# Definição do tamanho e origem da grade (ajuste conforme necessário)
grid_dim = [64, 64, 64]  # Grid dimensions (x, y, z)
grid_size = 1.0  # Size of each voxel in Ångström
origin = np.array([17.773, 63.285, 121.743])  # Origin of the grid in the PDB coordinate space


# Chamada de teste
parsed_pdb = parse_pdb("./3c9t.pdb")  # Substitua "your_pdb_file.pdb" pelo caminho do seu arquivo PDB
voxel_grid, info_grid = pdb_to_voxel_with_info(parsed_pdb)
plot_voxel(voxel_grid)



"""armazenar informações adicionais sobre cada átomo quando estiver preenchendo sua grade de voxels"""