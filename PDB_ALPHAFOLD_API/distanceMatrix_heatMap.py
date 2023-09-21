import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns

# Função para processar cada linha do arquivo CSV
def process_line_v5(line, total_columns):
    """Extract the sequence name and the distance values from a line using the provided rule."""
    parts = line.split(",")
    sequence_name = parts[0]
    values = []
    i = 1
    while i < len(parts):
        # Se o valor atual é 0 e é seguido por uma vírgula, é um valor zero completo
        if parts[i] == '0' and i < len(parts) - 1 and parts[i+1] == '0':
            values.append(0.0)
            i += 1
        # Se o valor atual é 0 e é seguido por um número, é um valor decimal
        elif parts[i] == '0' and i < len(parts) - 1:
            values.append(float(f"0.{parts[i+1]}"))
            i += 2
        else:
            i += 1
            
    # Garante que o número de valores corresponda ao total de colunas
    while len(values) < total_columns:
        values.append(0.0)

    return sequence_name, values

# Determine o número de colunas pelo número de nomes de sequências
total_columns = len(sequence_names_v4)

sequence_names_v5 = []
data_v5 = []

with open('/distance_matrix.csv', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        sequence_name, values = process_line_v5(line, total_columns)
        sequence_names_v5.append(sequence_name)
        data_v5.append(values)

# Construa o DataFrame final usando os novos dados
final_distance_matrix_v6 = pd.DataFrame(data_v5, columns=sequence_names_v5, index=sequence_names_v5)

# Remova a primeira linha que contém apenas zeros
final_distance_matrix_v7 = final_distance_matrix_v6.drop(final_distance_matrix_v6.index[0])

# Visualize a matriz de distância atualizada como um heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(final_distance_matrix_v7.iloc[:11, :11], cmap="YlGnBu", annot=True, fmt=".3f", linewidths=.5)
plt.title("Matriz de Distância Atualizada (Primeiros 11x11 Valores)")
plt.show()
