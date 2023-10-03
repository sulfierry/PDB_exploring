import pandas as pd

def read_csv_in_chunks(filename, chunk_size=50000):
    """Lê um arquivo CSV em chunks."""
    chunks = pd.read_csv(filename, sep=';', chunksize=chunk_size)
    return chunks

def filter_kinase_targets(chunk):
    """Filtrar moléculas que têm proteínas quinases como alvo."""
    kinase_chunk = chunk[chunk['Targets'].str.contains('kinase', case=False, na=False)]
    return kinase_chunk

def get_unique_kinases(kinase_chunk):
    """Obter o nome de todas as quinases que são alvos."""
    kinases = kinase_chunk['Targets'].unique()
    return kinases

def count_molecules_per_kinase(kinase_chunk):
    """Contar quantas moléculas por alvo de quinase existem."""
    kinase_counts = kinase_chunk['Targets'].value_counts()
    return kinase_counts

def main():
    filename = 'path_to_your_file.csv'
    chunks = read_csv_in_chunks(filename)
    
    all_kinases = []
    all_counts = {}

    for chunk in chunks:
        kinase_chunk = filter_kinase_targets(chunk)
        unique_kinases = get_unique_kinases(kinase_chunk)
        all_kinases.extend(unique_kinases)
        counts = count_molecules_per_kinase(kinase_chunk)
        for kinase, count in counts.items():
            if kinase in all_counts:
                all_counts[kinase] += count
            else:
                all_counts[kinase] = count

    # Imprimindo resultados
    print("Unique kinases:", list(set(all_kinases)))
    print("\nNumber of molecules per kinase target:")
    for kinase, count in all_counts.items():
        print(kinase, ":", count)

if __name__ == "__main__":
    main()
