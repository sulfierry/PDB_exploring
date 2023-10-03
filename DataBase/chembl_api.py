from chembl_webresource_client.new_client import new_client

def fetch_kinase_related_compounds():
    # Inicialize o cliente para a entidade 'molecule'
    molecule = new_client.molecule
    
    # Defina os termos de busca relacionados a inibidores de cinase
    search_terms = [
        "kinase inhibitor",
        "tyrosine kinase inhibitor",
        "TKI",
        "AKI",
        "MKI",
        "kinase antagonist",
        "EGFR kinase inhibitor",
    ]
    
    # Busque moléculas usando cada termo e coleciona resultados únicos
    unique_molecules = set()
    for term in search_terms:
        kinase_related = molecule.filter(description__icontains=term)
        for compound in kinase_related:
            unique_id = compound["molecule_chembl_id"]
            if unique_id not in unique_molecules:
                unique_molecules.add(unique_id)
                yield {
                    "chembl_id": unique_id,
                    "name": compound.get("pref_name", None) or (compound["molecule_synonyms"][0]["synonyms"] if compound["molecule_synonyms"] else "N/A"),
                    "description": compound.get("description", "N/A"),
                    "type": compound.get("molecule_type", "N/A")
                }

if __name__ == "__main__":
    kinase_related_drugs = list(fetch_kinase_related_compounds())
    for drug in kinase_related_drugs:
        print(f"ChEMBL ID: {drug['chembl_id']}, Name: {drug['name']}, Type: {drug['type']}, Description: {drug['description']}")
