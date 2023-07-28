import sys
from Bio.PDB import PDBParser, PDBIO

# Get the input and output file names from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Create a parser
parser = PDBParser(QUIET=True)

# Read the PDB file
structure = parser.get_structure("name", input_file)

# Create a PDBIO object and set the structure to be its model
io = PDBIO()
io.set_structure(structure)

# Save the structure to a new PDB file
io.save(output_file)
