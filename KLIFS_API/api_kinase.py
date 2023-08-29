import requests
import csv

# Lista de URLs
urls = [
    'https://klifs.net/api_v2/kinase_groups',         # groups
    'https://klifs.net/api_v2/kinase_families',       # families
    'https://klifs.net/api_v2/kinase_names',          # kinase names
    'https://klifs.net/api_v2/kinase_information',    # kinase information
    'https://klifs.net/api_v2/ligands_list',          # ligand list
    'https://klifs.net/api_v2/drug_list',             # drug list
    'https://klifs.net/api_v2/interactions_get_types' # interaction get types

    
]

def save_to_csv(data, headers, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for item in data:
            writer.writerow(item)

def print_and_save_dict_format(data, output_name):
    headers = list(data[0].keys())
    column_widths = [max(len(str(item[header])) for item in data + [dict.fromkeys(headers, header)]) for header in headers]
    
    header_row = "  ".join([headers[i].ljust(column_widths[i]) for i in range(len(headers))])
    print(header_row)
    print("-" * len(header_row))
    
    for item in data:
        row = [str(item[header]).ljust(column_widths[i]) for i, header in enumerate(headers)]
        print("  ".join(row))
    
    save_to_csv(data, headers, f'{output_name}.csv')

def print_and_save_list_format(data, output_name):
    for item in data:
        print(item)
    
    with open(f'{output_name}.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for item in data:
            writer.writerow([item])

for url in urls:
    # Extrair a última palavra da URL após "/"
    output_name = url.rstrip('/').split('/')[-1]

    response = requests.get(url, headers={'accept': 'application/json'})

    if response.status_code == 200:
        data = response.json()
        
        if isinstance(data, list):
            if isinstance(data[0], dict):
                print_and_save_dict_format(data, output_name)
            else:
                print_and_save_list_format(data, output_name)
        else:
            print("Formato de dados desconhecido.")
        
        print(f"\nDados salvos em '{output_name}.csv'")
        
    else:
        print(f"Erro ao acessar a API. Código de status: {response.status_code}")
        print("Mensagem de erro:", response.text)


