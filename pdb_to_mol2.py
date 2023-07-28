#python3.9 script.py input.pdb output.mol2

from rdkit import Chem
import sys


mol = Chem.MolFromPDBFile(sys.argv[1], removeHs=False)
Chem.MolToMolFile(mol, sys.argv[2] )
