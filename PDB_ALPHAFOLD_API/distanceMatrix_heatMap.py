# Consolidating the code into one block

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Function to process each line of the CSV file
def process_line_v5(line, total_columns):
    """Extract the sequence name and the distance values from a line using the provided rule."""
    parts = line.split(",")
    sequence_name = parts[0]
    values = []
    i = 1
    while i < len(parts):
        if parts[i] == '0' and i < len(parts) - 1 and parts[i+1] == '0':
            values.append(0.0)
            i += 1
        elif parts[i] == '0' and i < len(parts) - 1:
            values.append(float(f"0.{parts[i+1]}"))
            i += 2
        else:
            i += 1
    while len(values) < total_columns:
        values.append(0.0)
    return sequence_name, values

# Read the CSV file to determine the number of columns
with open('./distance_matrix.csv', 'r') as f:
    sequence_names = [line.strip().split(",")[0] for line in f if line.strip()]

total_columns = len(sequence_names)

data = []

with open('./distance_matrix.csv', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        sequence_name, values = process_line_v5(line, total_columns)
        data.append(values)

# Construct the DataFrame using the processed data
distance_matrix = pd.DataFrame(data, columns=sequence_names, index=sequence_names)

# Remove the first row which contains only zeros
distance_matrix = distance_matrix.drop(distance_matrix.index[0])

# Display the first 10x10 values of the matrix as a heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(distance_matrix.iloc[:11, :11], cmap="YlGnBu", annot=True, fmt=".3f", linewidths=.5)
plt.title("Matriz de Distância (Primeiros 10x10 Valores)")
plt.show()

