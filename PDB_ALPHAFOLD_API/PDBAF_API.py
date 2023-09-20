import os
import csv
import requests

class PDB_AF:
    ALPHAFOLD_API_ENDPOINT = "https://www.alphafold.ebi.ac.uk/api/prediction/"

    def __init__(self, pdb_id=None, uniprot_id=None):
        self.pdb_id = pdb_id
        self.uniprot_id = uniprot_id

    def download_pdb(self):
        url = f"https://files.rcsb.org/download/{self.pdb_id}.pdb"
        response = requests.get(url)

        if response.status_code == 200:
            with open(f"{self.pdb_id}.pdb", 'w') as file:
                file.write(response.text)
            print(f"Arquivo {self.pdb_id}.pdb baixado com sucesso!")
        else:
            print(f"Não foi possível baixar o PDB com o ID: {self.pdb_id}. Verifique se o ID é válido.")

    def fetch_alfafold_data(self):
        response = requests.get(self.ALPHAFOLD_API_ENDPOINT + self.uniprot_id)

        if response.status_code != 200:
            print(f"Erro ao acessar a API para o UniProt ID: {self.uniprot_id}. Código de status: {response.status_code}")
            return None

        return response.json()

    def download_alfafold_pdb(self):
        data = self.fetch_alfafold_data()

        if data:
            pdb_url = data[0].get("pdbUrl", "")
            if pdb_url:
                self.download_pdb_file(pdb_url, f"{self.uniprot_id}.pdb")
            else:
                print("URL do arquivo PDB não encontrado para este UniProt ID.")

    def download_pdb_file(self, pdb_url, filename):
        if os.path.exists(filename):
            print(f"O arquivo {filename} já foi baixado. Como os PDBs são iguais, somente uma cópia será salva.")
            return

        response = requests.get(pdb_url)
        
        if response.status_code == 200:
            with open(filename, 'w') as file:
                file.write(response.text)
            print(f"Arquivo PDB baixado com sucesso como {filename}!")
        else:
            print(f"Erro ao baixar o arquivo PDB. Código de status: {response.status_code}")



    @classmethod
    def process_csv(cls, filename):
        with open(filename, 'r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header row
            for row in reader:
                api, structure_id = row
                if api == "pdb":
                    pdb_instance = cls(pdb_id=structure_id)
                    pdb_instance.download_pdb()
                elif api == "af":
                    af_instance = cls(uniprot_id=structure_id)
                    af_instance.download_alfafold_pdb()
                else:
                    print(f"API desconhecida: {api}")

if __name__ == "__main__":
     

     PDB_AF.process_csv("example.csv")

    # pdb = PDB_AF(pdb_id="3FD5", uniprot_id="O15067")
    # pdb.download_pdb()
    # pdb.download_alfafold_pdb()
