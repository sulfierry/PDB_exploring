import requests

class PDBAF_API:
    ALPHAFOLD_API_ENDPOINT = "https://www.alphafold.ebi.ac.uk/api/prediction/"

    @staticmethod
    def pdb_API(pdb_id):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)

        if response.status_code == 200:
            with open(f"{pdb_id}.pdb", 'w') as file:
                file.write(response.text)
            print(f"Arquivo {pdb_id}.pdb baixado com sucesso!")
        else:
            print(f"Não foi possível baixar o PDB com o ID: {pdb_id}. Verifique se o ID é válido.")

    @staticmethod
    def af_API(uniprot_accession):
        response = requests.get(PDBAF_API.ALPHAFOLD_API_ENDPOINT + uniprot_accession)

        if response.status_code != 200:
            print(f"Erro ao acessar a API para o UniProt ID: {uniprot_accession}. Código de status: {response.status_code}")
            return None

        return response.json()

    @staticmethod
    def download_pdb_file(pdb_url, filename):
        response = requests.get(pdb_url)
        
        if response.status_code == 200:
            with open(filename, 'w') as file:
                file.write(response.text)
            print(f"Arquivo PDB baixado com sucesso como {filename}!")
        else:
            print(f"Erro ao baixar o arquivo PDB. Código de status: {response.status_code}")


if __name__ == "__main__":
    # Instancia da classe não é necessária, pois estamos usando métodos estáticos.
    
    pdb_id = "3FD5"
    PDBAF_API.pdb_API(pdb_id)

    uniprot_id = "O15067"
    data = PDBAF_API.af_API(uniprot_id)

    if data:
        pdb_url = data[0].get("pdbUrl", "")
        if pdb_url:
            PDBAF_API.download_pdb_file(pdb_url, f"{uniprot_id}.pdb")
        else:
            print("URL do arquivo PDB não encontrado para este UniProt ID.")
