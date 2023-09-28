import requests

# Defina o CID obtido
cid = 11531745

# Construa a URL para buscar compostos similares no PubChem
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON"

# Realiza a solicitação GET para a API
response = requests.get(url)

# Verifique se a solicitação foi bem-sucedida
if response.status_code == 200:
    data = response.json()
    similar_cids = data['IdentifierList']['CID']
    
    # Imprima os CIDs das moléculas similares
    for similar_cid in similar_cids:
        print(similar_cid)
else:
    print(f"Erro na solicitação HTTP. Código de status: {response.status_code}")
