import sys

def format_pdb(pdb_file, output_file):

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Filter out lines we don't care about
    lines = [line for line in lines if line.startswith(("ATOM", "HETATM", "TER"))]

    # Prepare to format PDB
    hetatm_flag = False  # Flag to track HETATM entries
    atom_counter = 1
    residue_counter = 1
    prev_res_num = None

    formatted_lines = []

    for line in lines:
        # Check if the current line is a HETATM or after TER
        if line.startswith("TER"):
            hetatm_flag = True
            formatted_lines.append(line)
            continue

        # Extract residue number from the line
        res_num = int(line[22:26].strip())

        # If we encounter a new residue, increment the residue counter
        if prev_res_num is not None and res_num != prev_res_num:
            residue_counter += 1

        # Assign residue counter for HETATM
        if hetatm_flag:
            line = line[:22] + f"{residue_counter:4}" + line[26:]

        # Adjust the line type to HETATM after TER
        if hetatm_flag:
            line = "HETATM" + line[6:]

        # Format atom number and residue number
        line = line[:6] + f"{atom_counter:5}" + line[11:22] + f"{residue_counter:4}" + line[26:]

        # Append the formatted line
        formatted_lines.append(line)

        # Update the atom counter and set the previous residue number
        atom_counter += 1
        prev_res_num = res_num

    # Add TER and END to the end of the file
    formatted_lines.append("TER\n")
    formatted_lines.append("END\n")

    with open(output_file, 'w') as f:
        f.writelines(formatted_lines)

    return formatted_lines

# Example usage:
# format_pdb("input.pdb", "output.pdb")


format_pdb(sys.argv[1], sys.argv[2])
