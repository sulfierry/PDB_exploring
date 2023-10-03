import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def compute_tanimoto_similarity(smiles, query_fp):
    """
    Compute Tanimoto similarity between a SMILES and a query fingerprint.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprint(mol, 2)
        return DataStructs.TanimotoSimilarity(query_fp, fp)
    else:
        return 0

def filter_tanimoto_similar_molecules(df, query_fp, threshold=0.8):
    """
    Filter molecules in a DataFrame based on Tanimoto similarity to a query fingerprint.
    """
    similar_molecules = df[df['SMILES'].apply(lambda x: compute_tanimoto_similarity(x, query_fp) >= threshold)]
    return similar_molecules

# Load the datasets
ligands_df = pd.read_csv('./klifs_ligands_list.csv')
drug_df = pd.read_csv('./klifs_drug_list.csv')


from rdkit import Chem

simplified_smiles = "c1ccc(nc1)Nc2ncc(c(n2)c3cc4c(c(c3)F)nc(n4)C)F"
simplified_mol = Chem.MolFromSmiles(simplified_smiles)

if simplified_mol:
    print("O SMILES é válido!")
else:
    print("O SMILES não é válido.")


# Define the query SMILES
query_smiles = simplified_smiles
query_mol = Chem.MolFromSmiles(query_smiles)
query_fp = AllChem.GetMorganFingerprint(query_mol, 2)

# Filter the ligands and drug datasets based on Tanimoto similarity
tanimoto_similar_ligands = filter_tanimoto_similar_molecules(ligands_df, query_fp)
tanimoto_similar_drugs = filter_tanimoto_similar_molecules(drug_df, query_fp)

# Save the filtered results to new CSV files
tanimoto_similar_ligands.to_csv('./tanimoto_similar_ligands.csv', index=False)
tanimoto_similar_drugs.to_csv('./tanimoto_similar_drugs.csv', index=False)
