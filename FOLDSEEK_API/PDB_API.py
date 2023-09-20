import requests

def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    # Verifica se a resposta é válida (código de status 200 indica sucesso)
    if response.status_code == 200:
        with open(f"{pdb_id}.pdb", 'w') as file:
            file.write(response.text)
        print(f"Arquivo {pdb_id}.pdb baixado com sucesso!")
    else:
        print(f"Não foi possível baixar o PDB com o ID: {pdb_id}. Verifique se o ID é válido.")

if __name__ == "__main__":
    pdb_id = ("3FD5")
    download_pdb(pdb_id)
