import requests
import csv
from chembl_webresource_client.new_client import new_client

def search_chembl(smiles):
    molecule = new_client.molecule
    similar_molecules = molecule.filter(molecule_structures__canonical_smiles__flexmatch=smiles)
    results_chembl = []
    for mol in similar_molecules:
        chembl_id = mol['molecule_chembl_id']
        canonical_smiles = mol['molecule_structures']['canonical_smiles']
        results_chembl.append(('ChEMBL', chembl_id, canonical_smiles))
    return results_chembl

def search_pubchem(smiles):
    url_smiles_to_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    response = requests.get(url_smiles_to_cid)
    if response.status_code != 200:
        return []
    cid = response.json()['IdentifierList']['CID'][0]
    
    url_similar = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON"
    response = requests.get(url_similar)
    if response.status_code != 200:
        return []
    
    similar_cids = response.json()['IdentifierList']['CID']
    results_pubchem = []
    for similar_cid in similar_cids:
        url_cid_to_smiles = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{similar_cid}/property/CanonicalSMILES/JSON"
        response = requests.get(url_cid_to_smiles)
        if response.status_code == 200:
            canonical_smiles = response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
            results_pubchem.append(('PubChem', similar_cid, canonical_smiles))
        else:
            results_pubchem.append(('PubChem', similar_cid, ''))

    return results_pubchem

import xml.etree.ElementTree as ET

def search_bindingdb(smiles, similarity_cutoff=0.8):
    base_url = "https://bindingdb.org/axis2/services/BDBService/getTargetByCompound"
    params = {
        'smiles': smiles,
        'cutoff': similarity_cutoff
    }

    response = requests.get(base_url, params=params)
    print(f"Response status: {response.status_code}")  # Log the status code
    print(f"Response content: {response.content}")  # Log the content

    if response.status_code != 200:
        return []

    # O BindingDB retorna os resultados em formato XML, então usaremos a biblioteca xml.etree.ElementTree para analisar
    root = ET.fromstring(response.content)

    results_bindingdb = []
    for compound in root.findall(".//compound"):
        smiles = compound.find("smiles").text
        for target in compound.findall("target"):
            uniprot = target.find("uniprot").text
            affinity = target.find("affinity").text
            affinity_type = target.find("affinity_type").text
            results_bindingdb.append(("BindingDB", smiles, uniprot, affinity_type, affinity))

    return results_bindingdb

def search_klifs(smiles, similarity_cutoff=0.85):
    base_url = "https://bindingdb.org/axis2/services/BDBService/getTargetByCompound"
    params = {
        'smiles': smiles,
        'cutoff': similarity_cutoff
    }

    response = requests.get(base_url, params=params)

    if response.status_code != 200:
        return []

    # O KLIFS retorna os resultados em formato XML, então usaremos a biblioteca xml.etree.ElementTree para analisar
    root = ET.fromstring(response.content)

    results_klifs = []
    for compound in root.findall(".//compound"):
        smiles = compound.find("smiles").text
        for target in compound.findall("target"):
            uniprot = target.find("uniprot").text
            affinity = target.find("affinity").text
            affinity_type = target.find("affinity_type").text
            results_klifs.append(("KLIFS", smiles, uniprot, affinity_type, affinity))

    return results_klifs

def save_to_csv(data, filename):
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Database', 'SMILES', 'ID/UniProt', 'Affinity Type', 'Affinity'])  # header
        for row in data:
            csvwriter.writerow(row)


def search_all_databases(smiles):
    results = []
    # results.extend(search_chembl(smiles))
    # results.extend(search_pubchem(smiles))
    results.extend(search_bindingdb(smiles))
    # results.extend(search_klifs(smiles))
    return results

def save_to_csv(results, filename='output.csv'):
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['NomeDB', 'CodigoDB', 'CodigoSmiles'])  # header
        for row in results:
            csvwriter.writerow(row)
            print(row)

def save_to_csv2(data, filename):
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Database', 'SMILES', 'ID/UniProt', 'Affinity Type', 'Affinity'])  # header
        for row in data:
            csvwriter.writerow(row)


query_smiles = "N[C@@H](CC(=O)Nc1ccc(Oc2cc(F)c(F)cc2Br)cc1)C(=O)O"
all_results = search_all_databases(query_smiles)
save_to_csv(all_results)


# results_klifs = search_klifs(query_smiles)
# save_to_csv2(results_klifs, 'KLIFS.csv')