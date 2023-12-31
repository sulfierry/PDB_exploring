from openbabel import openbabel
import sys

def convert_mol2_to_pdb(input_mol2_file, output_pdb_file):
    # Cria um objeto Open Babel para ler o arquivo .mol2
    ob_conversion = openbabel.OBConversion()
    ob_mol = openbabel.OBMol()

    # Carrega o arquivo .mol2
    if not ob_conversion.ReadFile(ob_mol, input_mol2_file):
        print(f"Erro ao ler o arquivo {input_mol2_file}")
        return

    # Cria um objeto Open Babel para escrever o arquivo .pdb
    ob_conversion.SetOutFormat("pdb")
    ob_conversion.WriteFile(ob_mol, output_pdb_file)

    print("Conversão concluída com sucesso.")

# Verifica se os argumentos foram fornecidos corretamente
if len(sys.argv) != 3:
    print("Uso: python script.py arquivo_input.mol2 arquivo_output.pdb")
else:
    input_mol2_file = sys.argv[1]
    output_pdb_file = sys.argv[2]

    convert_mol2_to_pdb(input_mol2_file, output_pdb_file)




def map_to_molecular_group(atom_name, residue_name):
    # Simplified mapping logic based on atom name, this needs refinement based on specific rules
    if atom_name.startswith("C"):
        return "CR"
    elif atom_name.startswith("O"):
        return "OR"
    elif atom_name.startswith("N"):
        return "NR"
    elif atom_name.startswith("H"):
        return "HR"
    elif atom_name.startswith("P"):
        return "PR"
    elif atom_name.startswith("F"):
        return "FR"

    return None
