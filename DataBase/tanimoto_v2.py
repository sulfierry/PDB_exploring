import pandas as pd
import re
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import csv

def find_smiles_column(df):
    """Identify the column which contains 'smiles' pattern in its name."""
    pattern = re.compile(r'.*smiles.*', re.IGNORECASE)
    for col in df.columns:
        if pattern.match(col):
            return col
    raise ValueError("No SMILES column found!")

def compute_tanimoto_similarity(smiles, query_fp):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return DataStructs.TanimotoSimilarity(query_fp, fp)
    else:
        return None

def filter_tanimoto_similar_molecules(df, query_fp, threshold=0.8):
    smiles_column = find_smiles_column(df)
    similar_molecules = df[df[smiles_column].apply(lambda x: (compute_tanimoto_similarity(x, query_fp) or 0) >= threshold)]
    return similar_molecules

# Load the datasets
klifs_ligands_df = pd.read_csv('./WAY_similar/KLIFS/klifs_ligands_list.csv')
klifs_drug_df = pd.read_csv('./WAY_similar/KLIFS/klifs_drug_list.csv')
pkid_drug_df = pd.read_csv('./WAY_similar/PKIDB/pkidb_2023-06-30.tsv', sep='\t')



simplified_smiles = "C1=CC(=CC=C1NC(=O)CC(C(=O)O)N)OC2=CC(=C(C=C2Br)F)F"
simplified_mol = Chem.MolFromSmiles(simplified_smiles)

if not simplified_mol:
    raise ValueError("O SMILES não é válido.")

query_mol = simplified_mol
query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)

tanimoto_similar_ligands = filter_tanimoto_similar_molecules(klifs_ligands_df, query_fp)
tanimoto_similar_drugs = filter_tanimoto_similar_molecules(klifs_drug_df, query_fp) 
tanimoto_similar_ligands = filter_tanimoto_similar_molecules(pkid_drug_df, query_fp)

tanimoto_similar_ligands.to_csv('klifs_tanimoto_similar_ligands.csv', index=False)
tanimoto_similar_drugs.to_csv('klifs_tanimoto_similar_drugs.csv', index=False)
tanimoto_similar_drugs.to_csv('pkid_tanimoto_similar_drugs.csv', index=False)
