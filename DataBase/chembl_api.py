from chembl_webresource_client.new_client import new_client

def fetch_all_compounds():
    molecule = new_client.molecule
    all_molecules = []

    # Definindo os parâmetros de paginação
    page_size = 100
    skip = 0

    while True:
        # Obtendo compostos em lotes
        batch = molecule.filter(molecule_type="Small molecule")[skip:skip + page_size]
        if not batch:
            break

        all_molecules.extend(batch)
        skip += page_size

    return all_molecules

def main():
    compounds = fetch_all_compounds()

    # Salvando os dados em um arquivo JSON
    with open('chembl_data.json', 'w') as f:
        import json
        json.dump(compounds, f)

if __name__ == "__main__":
    main()
