#                Copyright (c) Merck and Co., Inc., 1994
#                         All Rights Reserved

# Transcribe by: Leon Sulfierry GitHub:https://github.com/sulfierry

molecular_group_dict = {

    'CR': {'PRIMARY MMF TYPE': '1', 'DEFAULT TYPES': ['1', '1', '1', '0'], 'DEFINITION': 'ALKYL CARBON'},
    'C=C': {'PRIMARY MMF TYPE': '2', 'DEFAULT TYPES': ['2', '2', '1', '0'], 'DEFINITION': 'VINYLIC'},
    'C=O': {'PRIMARY MMF TYPE': '3', 'DEFAULT TYPES': ['3', '3', '1', '0'], 'DEFINITION': 'GENERAL CARBONYL C'},
    'CSP': {'PRIMARY MMF TYPE': '4', 'DEFAULT TYPES': ['4', '4', '1', '0'], 'DEFINITION': 'ACETYLENIC C'},
    'HC': {'PRIMARY MMF TYPE': '5', 'DEFAULT TYPES': ['5', '5', '5', '0'], 'DEFINITION': 'H-C'},
    'OR': {'PRIMARY MMF TYPE': '6', 'DEFAULT TYPES': ['6', '6', '6', '0'], 'DEFINITION': 'O-CSP3'},
    'O=C': {'PRIMARY MMF TYPE': '7', 'DEFAULT TYPES': ['7', '7', '6', '0'], 'DEFINITION': 'O=C, GENERIC'},
    'NR': {'PRIMARY MMF TYPE': '8', 'DEFAULT TYPES': ['8', '8', '8', '0'], 'DEFINITION': 'AMINE N'},
    'N=C': {'PRIMARY MMF TYPE': '9', 'DEFAULT TYPES': ['9', '9', '8', '0'], 'DEFINITION': 'N=C, IMINES'},
    'NC=O': {'PRIMARY MMF TYPE': '10', 'DEFAULT TYPES': ['10', '10', '8', '0'], 'DEFINITION': 'N-C=O, AMIDES'},
    'F': {'PRIMARY MMF TYPE': '11', 'DEFAULT TYPES': ['11', '11', '11', '0'], 'DEFINITION': 'FLUORINE'},
    'CL': {'PRIMARY MMF TYPE': '12', 'DEFAULT TYPES': ['12', '12', '12', '0'], 'DEFINITION': 'CHLORINE'},
    'BR': {'PRIMARY MMF TYPE': '13', 'DEFAULT TYPES': ['13', '13', '13', '0'], 'DEFINITION': 'BROMINE'},
    'I': {'PRIMARY MMF TYPE': '14', 'DEFAULT TYPES': ['14', '14', '14', '0'], 'DEFINITION': 'IODINE'},
    'S': {'PRIMARY MMF TYPE': '15', 'DEFAULT TYPES': ['15', '15', '15', '0'], 'DEFINITION': 'THIOL, SULFIDE'},
    'S=C': {'PRIMARY MMF TYPE': '16', 'DEFAULT TYPES': ['16', '16', '15', '0'], 'DEFINITION': 'S DOUBLY BONDED TO C'},
    'S=O': {'PRIMARY MMF TYPE': '17', 'DEFAULT TYPES': ['17', '17', '15', '0'], 'DEFINITION': 'SULFOXIDE S'},
    'SO2': {'PRIMARY MMF TYPE': '18', 'DEFAULT TYPES': ['18', '18', '15', '0'], 'DEFINITION': 'SULFONE S'},
    'SI': {'PRIMARY MMF TYPE': '19', 'DEFAULT TYPES': ['19', '19', '19', '0'], 'DEFINITION': 'SILICON'},
    'CR4R': {'PRIMARY MMF TYPE': '20', 'DEFAULT TYPES': ['20', '1', '1', '0'], 'DEFINITION': 'C IN CYCLOBUTYL'},
    'HOR': {'PRIMARY MMF TYPE': '21', 'DEFAULT TYPES': ['21', '21', '5', '0'], 'DEFINITION': 'H-O, ALCOHOLS'},
    'CR3R': {'PRIMARY MMF TYPE': '22', 'DEFAULT TYPES': ['22', '22', '1', '0'], 'DEFINITION': 'C IN CYCLOPROPLY'},
    'HNR': {'PRIMARY MMF TYPE': '23', 'DEFAULT TYPES': ['23', '23', '5', '0'], 'DEFINITION': 'H-N, AMINES'},
    'HOCO': {'PRIMARY MMF TYPE': '24', 'DEFAULT TYPES': ['24', '24', '5', '0'], 'DEFINITION': 'H-O, ACIDS'},
    'PO4': {'PRIMARY MMF TYPE': '25', 'DEFAULT TYPES': ['25', '25', '25', '0'], 'DEFINITION': 'PHOSPHODIESTER'},
    'P': {'PRIMARY MMF TYPE': '26', 'DEFAULT TYPES': ['26', '26', '25', '0'], 'DEFINITION': 'TRICOORDINATE P'},
    'HN=C': {'PRIMARY MMF TYPE': '27', 'DEFAULT TYPES': ['27', '28', '5', '0'], 'DEFINITION': 'IMINE N-H'},
    'HNCO': {'PRIMARY MMF TYPE': '28', 'DEFAULT TYPES': ['28', '28', '5', '0'], 'DEFINITION': 'H-N, AMIDES'},
    'HOCC': {'PRIMARY MMF TYPE': '29', 'DEFAULT TYPES': ['29', '29', '5', '0'], 'DEFINITION': 'H-O, ENOLS, PHENOLS'},
    'CE4R': {'PRIMARY MMF TYPE': '30', 'DEFAULT TYPES': ['30', '2', '1', '0'], 'DEFINITION': 'C=C IN 4-RING'},
    'HOH': {'PRIMARY MMF TYPE': '31', 'DEFAULT TYPES': ['31', '31', '31', '0'], 'DEFINITION': 'H-OH'},
    'O2CM': {'PRIMARY MMF TYPE': '32', 'DEFAULT TYPES': ['32', '7', '6', '0'], 'DEFINITION': 'O, CARBOXYLATE ANION'},
    'HOS': {'PRIMARY MMF TYPE': '33', 'DEFAULT TYPES': ['33', '21', '5', '0'], 'DEFINITION': 'H-O-S, SULF ACIDS'},
    'NR+': {'PRIMARY MMF TYPE': '34', 'DEFAULT TYPES': ['34', '8', '8', '0'], 'DEFINITION': 'N+, QUATERNARY N'},
    'OM': {'PRIMARY MMF TYPE': '35', 'DEFAULT TYPES': ['35', '6', '6', '0'], 'DEFINITION': 'OXIDE OXYGEN ON SP3 C'},
    'HNR+': {'PRIMARY MMF TYPE': '36', 'DEFAULT TYPES': ['36', '36', '5', '0'], 'DEFINITION': 'H-N+'},
    'CB': {'PRIMARY MMF TYPE': '37', 'DEFAULT TYPES': ['37', '2', '1', '0'], 'DEFINITION': 'AROMATIC C'},
    'NPYD': {'PRIMARY MMF TYPE': '38', 'DEFAULT TYPES': ['38', '9', '8', '0'], 'DEFINITION': 'AROMATIC N, PYRIDINE'},
    'NPYL': {'PRIMARY MMF TYPE': '39', 'DEFAULT TYPES': ['39', '10', '8', '0'], 'DEFINITION': 'AROMATIC N, PYRROLE'},
    'NC=C': {'PRIMARY MMF TYPE': '40', 'DEFAULT TYPES': ['40', '10', '8', '0'], 'DEFINITION': 'N-C=C (DELOC LP)'},
    'CO2M': {'PRIMARY MMF TYPE': '41', 'DEFAULT TYPES': ['41', '3', '1', '0'], 'DEFINITION': 'C IN CO2- ANION'},
    'NSP': {'PRIMARY MMF TYPE': '42', 'DEFAULT TYPES': ['42', '42', '8', '0'], 'DEFINITION': 'N TRIPLE BONDED'},
    'NSO2': {'PRIMARY MMF TYPE': '43', 'DEFAULT TYPES': ['43', '10', '8', '0'], 'DEFINITION': 'N, SULFONAMIDES'},
    'STHI': {'PRIMARY MMF TYPE': '44', 'DEFAULT TYPES': ['44', '16', '15', '0'], 'DEFINITION': 'S IN THIOPHENE'},
    'NO2': {'PRIMARY MMF TYPE': '45', 'DEFAULT TYPES': ['45', '10', '8', '0'], 'DEFINITION': 'NITRO GROUP N'},
    'N=O': {'PRIMARY MMF TYPE': '46', 'DEFAULT TYPES': ['46', '9', '8', '0'], 'DEFINITION': 'NITROSO GROUP N'},
    'NAZT': {'PRIMARY MMF TYPE': '47', 'DEFAULT TYPES': ['47', '42', '8', '0'], 'DEFINITION': 'TERMINAL N, AZIDE'},
    'NSO': {'PRIMARY MMF TYPE': '48', 'DEFAULT TYPES': ['48', '9', '8', '0'], 'DEFINITION': 'DIVAL. N IN S(N)(O) GP'},
    'O+': {'PRIMARY MMF TYPE': '49', 'DEFAULT TYPES': ['49', '6', '6', '0'], 'DEFINITION': 'OXONIUM (TRICOORD) O'},
    'HO+': {'PRIMARY MMF TYPE': '50', 'DEFAULT TYPES': ['50', '21', '5', '0'], 'DEFINITION': 'H ON OXONIUM OXYGEN'},
    'O=+': {'PRIMARY MMF TYPE': '51', 'DEFAULT TYPES': ['51', '7', '6', '0'], 'DEFINITION': 'OXENIUM OXYGEN+'},
    'HO=+': {'PRIMARY MMF TYPE': '52', 'DEFAULT TYPES': ['52', '21', '5', '0'], 'DEFINITION': 'H ON OXENIUM O+'},
    '=N=': {'PRIMARY MMF TYPE': '53', 'DEFAULT TYPES': ['53', '42', '8', '0'], 'DEFINITION': 'N TWICE DOUBLE BONDED'},
    'N+=C': {'PRIMARY MMF TYPE': '54', 'DEFAULT TYPES': ['54', '9', '8', '0'], 'DEFINITION': 'IMINIUM NITROGEN'},
    'NCN+': {'PRIMARY MMF TYPE': '55', 'DEFAULT TYPES': ['55', '10', '8', '0'], 'DEFINITION': 'N IN +N=C-N: ; Q=1/2'},
    'NGD+': {'PRIMARY MMF TYPE': '56', 'DEFAULT TYPES': ['56', '10', '8', '0'], 'DEFINITION': 'GUANIDINIUM N; Q=1/3'},
    'CGD+': {'PRIMARY MMF TYPE': '57', 'DEFAULT TYPES': ['57', '2', '1', '0'], 'DEFINITION': 'GUANIDINIUM CARBON'},
    'NPD+': {'PRIMARY MMF TYPE': '58', 'DEFAULT TYPES': ['58', '10', '8', '0'], 'DEFINITION': 'N PYRIDINIUM ION'},
    'OFUR': {'PRIMARY MMF TYPE': '59', 'DEFAULT TYPES': ['59', '6', '6', '0'], 'DEFINITION': 'AROMATIC O, FURAN'},
    'C%': {'PRIMARY MMF TYPE': '60', 'DEFAULT TYPES': ['60', '4', '1', '0'], 'DEFINITION': 'ISONITRILE CARBON'},
    'NR%': {'PRIMARY MMF TYPE': '61', 'DEFAULT TYPES': ['61', '42', '8', '0'], 'DEFINITION': 'ISONITRILE N'},
    'NM': {'PRIMARY MMF TYPE': '62', 'DEFAULT TYPES': ['62', '10', '8', '0'], 'DEFINITION': 'SULFONAMIDE N-'},
    'C5A': {'PRIMARY MMF TYPE': '63', 'DEFAULT TYPES': ['63', '2', '1', '0'], 'DEFINITION': 'ALPHA AROM 5-RING C'},
    'C5B': {'PRIMARY MMF TYPE': '64', 'DEFAULT TYPES': ['64', '2', '1', '0'], 'DEFINITION': 'BETA AROM 5-RING C'},
    'N5A': {'PRIMARY MMF TYPE': '65', 'DEFAULT TYPES': ['65', '9', '8', '0'], 'DEFINITION': 'ALPHA AROM 5-RING N'},
    'N5B': {'PRIMARY MMF TYPE': '66', 'DEFAULT TYPES': ['66', '9', '8', '0'], 'DEFINITION': 'ALPHA AROM 5-RING N'},
    'N2OX': {'PRIMARY MMF TYPE': '67', 'DEFAULT TYPES': ['67', '9', '8', '0'], 'DEFINITION': 'NITROGEN IN N-OXIDE'},
    'N3OX': {'PRIMARY MMF TYPE': '68', 'DEFAULT TYPES': ['68', '8', '8', '0'], 'DEFINITION': 'NITROGEN IN N-OXIDE'},
    'NPOX': {'PRIMARY MMF TYPE': '69', 'DEFAULT TYPES': ['69', '9', '8', '0'], 'DEFINITION': 'NITROGEN IN N-OXIDE'},
    'OH2': {'PRIMARY MMF TYPE': '70', 'DEFAULT TYPES': ['70', '70', '70', '70'], 'DEFINITION': 'OXYGEN IN WATER'},
    'HS': {'PRIMARY MMF TYPE': '71', 'DEFAULT TYPES': ['71', '5', '5', '0'], 'DEFINITION': 'H-S'},
    'S2CM': {'PRIMARY MMF TYPE': '72', 'DEFAULT TYPES': ['72', '16', '15', '0'], 'DEFINITION': 'THIOCARBOXYLATE S'},
    'SO2M': {'PRIMARY MMF TYPE': '73', 'DEFAULT TYPES': ['73', '18', '15', '0'], 'DEFINITION': 'SULFUR IN SULFINATE'},
    '=S=O': {'PRIMARY MMF TYPE': '74', 'DEFAULT TYPES': ['74', '17', '15', '0'], 'DEFINITION': 'SULFINYL SULFUR, C=S=O'},
    '-P=C': {'PRIMARY MMF TYPE': '75', 'DEFAULT TYPES': ['75', '26', '25', '0'], 'DEFINITION': 'P DOUBLY BONDED TO C'},
    'N5M': {'PRIMARY MMF TYPE': '76', 'DEFAULT TYPES': ['76', '9', '8', '0'], 'DEFINITION': 'NEG N IN TETRAZOLE AN'},
    'CLO4': {'PRIMARY MMF TYPE': '77', 'DEFAULT TYPES': ['77', '12', '12', '0'], 'DEFINITION': 'CHLORINE IN CLO4(-)'},
    'C5': {'PRIMARY MMF TYPE': '78', 'DEFAULT TYPES': ['78', '2', '1', '0'], 'DEFINITION': 'GENERAL AROM 5-RING C'},
    'N5': {'PRIMARY MMF TYPE': '79', 'DEFAULT TYPES': ['79', '9', '8', '0'], 'DEFINITION': 'GENERAL AROM 5-RING N'},
    'CIM+': {'PRIMARY MMF TYPE': '80', 'DEFAULT TYPES': ['80', '2', '1', '0'], 'DEFINITION': 'C IN N-C-N, IM+ ION'},
    'NIM+': {'PRIMARY MMF TYPE': '81', 'DEFAULT TYPES': ['81', '10', '8', '0'], 'DEFINITION': 'N IN N-C-N, IM+ ION'},
    'N5AX': {'PRIMARY MMF TYPE': '82', 'DEFAULT TYPES': ['82', '9', '8', '0'], 'DEFINITION': '5R NITROGEN IN N-OXIDE'},
    'FE+2': {'PRIMARY MMF TYPE': '87', 'DEFAULT TYPES': ['87', '87', '87', '87'], 'DEFINITION': 'IRON +2 CATION'},
    'FE+3': {'PRIMARY MMF TYPE': '88', 'DEFAULT TYPES': ['88', '88', '88', '88'], 'DEFINITION': 'IRON +3 CATION'},
    'F-': {'PRIMARY MMF TYPE': '89', 'DEFAULT TYPES': ['89', '89', '89', '89'], 'DEFINITION': 'FLUORIDE ANION'},
    'CL-': {'PRIMARY MMF TYPE': '90', 'DEFAULT TYPES': ['90', '90', '90', '90'], 'DEFINITION': 'CHLORIDE ANION'},
    'BR-': {'PRIMARY MMF TYPE': '91', 'DEFAULT TYPES': ['91', '91', '91', '91'], 'DEFINITION': 'BROMIDE ANION'},
    'LI+': {'PRIMARY MMF TYPE': '92', 'DEFAULT TYPES': ['92', '92', '92', '92'], 'DEFINITION': 'LITHIUM CATION'},
    'NA+': {'PRIMARY MMF TYPE': '93', 'DEFAULT TYPES': ['93', '93', '93', '93'], 'DEFINITION': 'SODIUM CATION'},
    'K+': {'PRIMARY MMF TYPE': '94', 'DEFAULT TYPES': ['94', '94', '94', '94'], 'DEFINITION': 'POTASSIUM CATION'},
    'ZN+2': {'PRIMARY MMF TYPE': '95', 'DEFAULT TYPES': ['95', '95', '95', '95'], 'DEFINITION': 'DIPOSITIVE ZINC CATION'},
    'CA+2': {'PRIMARY MMF TYPE': '96', 'DEFAULT TYPES': ['96', '96', '96', '96'], 'DEFINITION': 'DIPOSITIVE CALCIUM CATION'},
    'CU+1': {'PRIMARY MMF TYPE': '97', 'DEFAULT TYPES': ['97', '97', '97', '97'], 'DEFINITION': 'MONOPOSITIVE COPPER CATION'},
    'CU+2': {'PRIMARY MMF TYPE': '98', 'DEFAULT TYPES': ['98', '98', '98', '98'], 'DEFINITION': 'DIPOSITIVE COPPER CATION'},
    'MG+2': {'PRIMARY MMF TYPE': '99', 'DEFAULT TYPES': ['99', '99', '99', '99'], 'DEFINITION': 'DIPOSITIVE MAGNESIUM CATION'},
    
}
