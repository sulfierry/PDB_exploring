import requests

def submit_to_foldseek(file_path):
    # URL para a API
    url = "https://search.foldseek.com/api/ticket"

    # Definição dos parâmetros e arquivos a serem enviados
    databases = [
        'afdb50', 
        'afdb-swissprot', 
        'afdb-proteome', 
        'cath50', 
        'mgnify_esm30', 
        'pdb100', 
        'gmgcl_id'
    ]
    files = {'q': (file_path, open(file_path, 'rb'))}
    data = {
        'mode': '3diaa',
    }
    for db in databases:
        data['database[]'] = db

    # Realizando a requisição POST
    response = requests.post(url, files=files, data=data)

    # Verifica se a requisição foi bem-sucedida
    if response.status_code == 200:
        print("Arquivo enviado com sucesso!")
        return response.json()
    else:
        print("Erro ao enviar o arquivo. Status code:", response.status_code)
        return None

# Teste
file_path = "YOUR_PATH_TO_FILE.pdb"
result = submit_to_foldseek(file_path)
if result:
    print(result)
