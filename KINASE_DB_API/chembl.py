from chembl_webresource_client.new_client import new_client

# Iniciar cliente de molécula
molecule = new_client.molecule

# Definir a estrutura química SMILES de sua query
query_smiles = "N[C@@H](CC(=O)Nc1ccc(Oc2cc(F)c(F)cc2Br)cc1)C(=O)O"

# Buscar moléculas similares
similar_molecules = molecule.filter(molecule_structures__canonical_smiles__flexmatch=query_smiles)

# Iterar sobre as moléculas similares encontradas
for mol in similar_molecules:
    print(mol['molecule_chembl_id'], mol['molecule_structures']['canonical_smiles'])

