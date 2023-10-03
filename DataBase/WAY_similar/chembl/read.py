import pandas as pd

# Definir o tamanho do batch
chunk_size = 50000  # por exemplo, ler 50.000 linhas de cada vez
filename =  'chembl_database_031023.csv'  # substitua por seu caminho de arquivo

import pandas as pd

def read_csv_in_chunks(filename, chunk_size=50000):
    """Lê um arquivo CSV em chunks."""
    chunks = pd.read_csv(filename, sep=';', chunksize=chunk_size)
    return chunks 

if __name__ == "__main__":
    filename = 'o.csv'
    chunks = read_csv_in_chunks(filename)
    
    for chunk in chunks:
        print(chunk.head())  # Imprime as primeiras linhas de cada chunk para verificação
