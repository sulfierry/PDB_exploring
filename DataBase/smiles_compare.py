from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors

# Exemplo de moléculas (use a sua lista completa)
molecules = [
    "CNc1cc(OC)c(C(N)=O)cn1",
    # ... (adicionar todas as moléculas)
]

# Estrutura para armazenar propriedades
data = []

for mol in molecules:
    molecule = Chem.MolFromSmiles(mol)
    if molecule:  # Se o SMILES for válido
        MW = Descriptors.MolWt(molecule)
        ClogP = Crippen.MolLogP(molecule)
        HBD = Lipinski.NumHDonors(molecule)
        HBA = Lipinski.NumHAcceptors(molecule)
        TPSA = rdMolDescriptors.CalcTPSA(molecule)
        NRB = Lipinski.NumRotatableBonds(molecule)
        NAR = rdMolDescriptors.CalcNumAromaticRings(molecule)
        NCA = rdMolDescriptors.CalcNumChiralCenters(molecule)
        
        data.append((mol, MW, ClogP, HBD, HBA, TPSA, NRB, NAR, NCA))

# Aqui, você pode inserir os dados na sua tabela ou filtrar os que atendem às suas condições.
