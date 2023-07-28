# python script.py protein.pdb ligand.pdb output.pdb

import sys

def get_last_atom_index(pdb_lines):
    # Função para extrair o último índice da segunda coluna do arquivo PDB
    last_index = 0
    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_index = int(line[6:11])
            last_index = max(last_index, atom_index)
    return last_index

def merge_pdbs(protein_file, ligand_file, output_file):
    # Ler as informações dos arquivos PDB da proteína e do ligante
    with open(protein_file, 'r') as protein_f, open(ligand_file, 'r') as ligand_f:
        protein_lines = protein_f.readlines()
        ligand_lines = ligand_f.readlines()

    # Obter o índice do último átomo da proteína
    last_protein_atom_index = get_last_atom_index(protein_lines)

    # Criar um novo arquivo PDB combinando as informações da proteína e do ligante
    with open(output_file, 'w') as output_f:
        # Copiar as informações da proteína para o novo arquivo
        for line in protein_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_index = int(line[6:11])
                # Atualizar o índice do átomo adicionando o deslocamento do último índice da proteína
                atom_index += last_protein_atom_index
                output_f.write(f"{line[:6]}{atom_index:5}{line[11:]}")
            elif line.startswith('END'):
                # Trocar 'END' por 'TER'
                output_f.write('TER\n')
            else:
                output_f.write(line)

        # Adicionar as informações do ligante ao novo arquivo, atualizando os índices dos átomos
        atom_offset = last_protein_atom_index
        ligand_index_offset = last_protein_atom_index + 1  # Índice do ligante deve começar após o último índice da proteína
        for line in ligand_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_index = int(line[6:11]) + ligand_index_offset
                output_f.write(f"{line[:6]}{atom_index:5}{line[11:]}")
            elif line.startswith('CONECT'):
                # Corrigir os índices no comando 'CONECT' para corresponderem aos índices do ligante
                indices = line.split()[1:]
                corrected_indices = [str(int(index) + ligand_index_offset) for index in indices]
                output_f.write('CONECT ' + ' '.join(corrected_indices) + '\n')
            else:
                output_f.write(line)

if __name__ == "__main__":
    # Verificar se todos os argumentos foram fornecidos
    if len(sys.argv) != 4:
        print("Uso: python merge_ligand_protein.py protein.pdb ligand.pdb output.pdb")
        sys.exit(1)

    # Obter os nomes dos arquivos de entrada e saída a partir dos argumentos da linha de comando
    protein_file = sys.argv[1]
    ligand_file = sys.argv[2]
    output_file = sys.argv[3]

    # Chamar a função para mesclar os arquivos PDB
    merge_pdbs(protein_file, ligand_file, output_file)
