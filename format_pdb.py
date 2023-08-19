import sys

AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]

def format_pdb(pdb_file, output_file):
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: {pdb_file} not found.")
        return

    lines = [line for line in lines if line.startswith(("ATOM", "HETATM", "TER"))]

    hetatm_flag = False
    atom_counter = 1
    residue_counter = 1
    prev_res_num = None

    formatted_lines = []

    for i, line in enumerate(lines):
        if line.startswith("TER"):
            # Check the residue type of the next line (if it exists)
            if i + 1 < len(lines) and lines[i + 1][17:20] in AMINO_ACIDS:
                hetatm_flag = False
            else:
                hetatm_flag = True
            formatted_lines.append(line)
            continue

        try:
            res_num = int(line[22:26].strip())
        except ValueError:
            print(f"Error parsing residue number in line: {line.strip()}")
            continue

        if prev_res_num is not None and res_num != prev_res_num:
            residue_counter += 1

        if hetatm_flag:
            line = line[:22] + f"{residue_counter:4}" + line[26:]
            line = "HETATM" + line[6:]
        else:
            line = line[:22] + f"{residue_counter:4}" + line[26:]
            line = "ATOM  " + line[6:]

        line = line[:6] + f"{atom_counter:5}" + line[11:22] + f"{residue_counter:4}" + line[26:]
        formatted_lines.append(line)

        atom_counter += 1
        prev_res_num = res_num

    formatted_lines.append("TER\n")
    formatted_lines.append("END\n")

    with open(output_file, 'w') as f:
        f.writelines(formatted_lines)

    return formatted_lines

if __name__ == "__main__":


    format_pdb(sys.argv[1], sys.argv[2])