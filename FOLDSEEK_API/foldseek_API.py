import requests
import time
from tqdm import tqdm

class FoldseekAPI:
    
    def __init__(self):
        self.MMSEQS2_BASE_URL = "https://search.mmseqs.com/api"
        self.FOLDSEEK_BASE_URL = "https://search.foldseek.com/api"

    def check_ticket_status(self, ticket_id):
        url = f"{self.MMSEQS2_BASE_URL}/ticket/{ticket_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            return response.json()
        else:
            print("Erro ao consultar o status. Status code:", response.status_code)
            return None

    def download_results(self, ticket_id):
        url = f"{self.MMSEQS2_BASE_URL}/result/download/{ticket_id}"
        response = requests.get(url, stream=True)
        with open(f"results_{ticket_id}.pdb", 'wb') as fd:
            for chunk in response.iter_content(chunk_size=128):
                fd.write(chunk)
        print("Download concluído.")

    def submit_to_foldseek(self, file_path):
        url = f"{self.FOLDSEEK_BASE_URL}/ticket"
        databases = [
            'afdb50', 'afdb-swissprot', 'afdb-proteome', 'cath50', 
            'mgnify_esm30', 'pdb100', 'gmgcl_id'
        ]

        data = {'mode': '3diaa'}
        for db in databases:
            data.setdefault('database[]', []).append(db)

        with open(file_path, 'rb') as f:
            files = {'q': (file_path, f)}
            response = requests.post(url, files=files, data=data)

        if response.status_code == 200:
            print("Arquivo enviado com sucesso!")
            return response.json()
        else:
            print("Erro ao enviar o arquivo. Status code:", response.status_code)
            return None

    def wait_for_completion(self, ticket_id, polling_interval=60):
        for _ in tqdm(range(0, 100), desc="Aguardando processamento", unit="poll"):
            status = self.check_ticket_status(ticket_id)
            if status and status.get('status') == "COMPLETED":
                print("\nProcessamento concluído!")
                return True
            elif status and status.get('status') == "ERROR":
                print("\nHouve um erro no processamento.")
                return False
            time.sleep(polling_interval)
        return False

# Exemplo de uso da classe:

api = FoldseekAPI()

# Teste
file_path = "./3c9t.pdb"
result = api.submit_to_foldseek(file_path)
if result:
    print(result)

ticket_id = result['id']
print(ticket_id)

# Aguarda a conclusão
if api.wait_for_completion(ticket_id):
    api.download_results(ticket_id)
else:
    print("O trabalho não foi concluído com sucesso.")
