import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import concurrent.futures
from tqdm import tqdm

# Definir o tamanho do batch
chunk_size = 50000
filename = 'o.csv'

def read_csv_in_chunks(filename, chunk_size=50000):
    """Lê um arquivo CSV em chunks."""
    chunks = pd.read_csv(filename, sep=';', chunksize=chunk_size)
    return chunks

def find_exact_smiles(chunk, smiles_query):
    """Busca uma correspondência exata para a query smiles no chunk."""
    exact_match = chunk[chunk['Smiles'] == smiles_query]
    return exact_match

def clean_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:  # Se não conseguir converter o SMILES
        return None
    Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol)

def tanimoto_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if not mol1 or not mol2:  # Se qualquer molécula for None
        return 0

    try:
        fp1 = AllChem.GetMorganFingerprint(mol1, 2)
        fp2 = AllChem.GetMorganFingerprint(mol2, 2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except:
        return 0  # Retorne 0 se houver um erro

def find_similar_smiles(chunk, smiles_query, threshold=0.85):
    """Busca smiles semelhantes à query no chunk usando a similaridade de Tanimoto."""
    # Remove as linhas que têm NaN na coluna 'Smiles'
    chunk = chunk.dropna(subset=['Smiles']).copy()  # Aqui, estamos fazendo uma cópia

    chunk['Smiles'] = chunk['Smiles'].apply(clean_molecule)  # Limpe os SMILES
    chunk = chunk.dropna(subset=['Smiles'])
    
    # Calcula a similaridade de Tanimoto para cada molécula no chunk em relação à query
    chunk['Similarity'] = chunk['Smiles'].apply(lambda x: tanimoto_similarity(smiles_query, x))
    
    # Filtra as moléculas com similaridade acima do limiar
    similar_molecules = chunk[chunk['Similarity'] > threshold]
    
    return similar_molecules.sort_values(by='Similarity', ascending=False)

def process_chunk(chunk, smiles_query):
    exact_results = find_exact_smiles(chunk, smiles_query)
    similar_results = find_similar_smiles(chunk, smiles_query)
    
    return exact_results, similar_results


if __name__ == "__main__":
    filename = 'o.csv'
    chunks = read_csv_in_chunks(filename)
    
    smiles_query = "CCOc1ccc(/C=N/NC(=O)c2ccc(/N=N/N(C)C)cc2)cc1"  # Exemplo de query
    
    # Número de workers (normalmente definido pelo número de CPUs)
    num_workers = 7 # Altere isso conforme o número de CPUs disponíveis em seu sistema
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # O método map do executor processará cada chunk em paralelo
        for exact_results, similar_results in tqdm(executor.map(process_chunk, chunks, [smiles_query]*chunk_size), desc="Processando chunks", unit="chunk"):
            if not exact_results.empty:
                print("Buscando por correspondência exata:")
                print(exact_results)
                
            if not similar_results.empty:
                print("\nBuscando por smiles similares:")
                print(similar_results)
