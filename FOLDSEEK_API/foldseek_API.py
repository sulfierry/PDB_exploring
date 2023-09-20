import requests

BASE_URL = "https://search.mmseqs.com/api"

def check_ticket_status(ticket_id):
    # URL para verificar o status do trabalho
    url = f"{BASE_URL}/ticket/{ticket_id}"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print("Erro ao consultar o status. Status code:", response.status_code)
        return None

def get_results(ticket_id, entry):
    # URL para obter os resultados
    url = f"{BASE_URL}/result/{ticket_id}/{entry}"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print("Erro ao obter os resultados. Status code:", response.status_code)
        return None

# Usando as funções
ticket_id = 'YOUR_TICKET_ID'
status = check_ticket_status(ticket_id)
if status and status.get('status') == "COMPLETED":
    # Assumindo 'entry' como um parâmetro necessário. Modifique conforme necessário.
    entry = 'YOUR_ENTRY'
    results = get_results(ticket_id, entry)
    print(results)
else:
    print("O trabalho ainda não está concluído ou ocorreu um erro.")




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


def check_foldseek_status(ticket_id):
    # URL para verificar o status do trabalho
    url = f"https://search.foldseek.com/api/status/{ticket_id}"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print("Erro ao consultar o status. Status code:", response.status_code)
        return None

import requests

BASE_URL = "https://search.mmseqs.com/api"

def check_ticket_status(ticket_id):
    # URL para verificar o status do trabalho
    url = f"{BASE_URL}/ticket/{ticket_id}"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print("Erro ao consultar o status. Status code:", response.status_code)
        return None

def get_results(ticket_id, entry):
    # URL para obter os resultados
    url = f"{BASE_URL}/result/{ticket_id}/{entry}"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print("Erro ao obter os resultados. Status code:", response.status_code)
        return None



# Teste
file_path = "./3c9t.pdb"
result = submit_to_foldseek(file_path)
if result:
    print(result)

print(result['id'])

# Usando a função para consultar o status
#status = check_foldseek_status(result['id'])
#if status:
#    print(status)



# Usando as funções
# ticket_id = 'YOUR_TICKET_ID'
status = check_ticket_status(result['id'])
if status and status.get('status') == "COMPLETED":
    # Assumindo 'entry' como um parâmetro necessário. Modifique conforme necessário.
    entry = 'YOUR_ENTRY'
    results = get_results(result['id'], entry)
    print(results)
else:
    print("O trabalho ainda não está concluído ou ocorreu um erro.")
