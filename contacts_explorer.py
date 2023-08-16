# Check possible contacts between structure and ligand
# python script.py ref.pdb LIG PDB_LIG.csv CHAIN(optional)
# GitHub: github.com/sulfierry/

import csv
import sys


# Variaveis globais ################################################################################################

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

vdw_dict = {
    '1': {'alpha-i': '1.050', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CR', 'Origin': 'E94'},
    '2': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C=C', 'Origin': 'E94'},
    '3': {'alpha-i': '1.100', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C=O', 'Origin': 'E94'},
    '4': {'alpha-i': '1.300', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CSP', 'Origin': 'E94'},
    '5': {'alpha-i': '0.250', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': '-', 'Symb': 'HC', 'Origin': 'C94'},
    '6': {'alpha-i': '0.70', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OR', 'Origin': 'C94'},
    '7': {'alpha-i': '0.65', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'O=C', 'Origin': 'C94'},
    '8': {'alpha-i': '1.15', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NR', 'Origin': 'C94'},
    '9': {'alpha-i': '0.90', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N=C', 'Origin': 'C94'},
    '10': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NC=O', 'Origin': 'E94'},
    '11': {'alpha-i': '0.35', 'N-i': '3.480', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'F', 'Origin': 'C94'},
    '12': {'alpha-i': '2.300', 'N-i': '5.100', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'CL', 'Origin': 'E94'},
    '13': {'alpha-i': '3.400', 'N-i': '6.000', 'A-i': '3.190', 'G-i': '1.359', 'DA': 'A', 'Symb': 'BR', 'Origin': 'E94'},
    '14': {'alpha-i': '5.500', 'N-i': '6.950', 'A-i': '3.080', 'G-i': '1.404', 'DA': 'A', 'Symb': 'I', 'Origin': 'E94'},
    '15': {'alpha-i': '3.00', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'S', 'Origin': 'C94'},
    '16': {'alpha-i': '3.900', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'S=C', 'Origin': 'E94'},
    '17': {'alpha-i': '2.700', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SO', 'Origin': 'E94'},
    '18': {'alpha-i': '2.100', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SO2', 'Origin': 'E94'},
    '19': {'alpha-i': '4.500', 'N-i': '4.200', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SI', 'Origin': 'E94'},
    '20': {'alpha-i': '1.050', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CR3R', 'Origin': 'E94'},
    '21': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOR', 'Origin': 'C94'},
    '22': {'alpha-i': '1.100', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CR3R', 'Origin': 'E94'},
    '23': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HNR', 'Origin': 'C94'},
    '24': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOCO', 'Origin': 'C94'},
    '25': {'alpha-i': '1.600', 'N-i': '4.500', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'PO4', 'Origin': 'E94'},
    '26': {'alpha-i': '3.600', 'N-i': '4.500', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'P', 'Origin': 'E94'},
    '27': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HN=C', 'Origin': 'C94'},
    '28': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HNCO', 'Origin': 'C94'},
    '29': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOCC', 'Origin': 'C94'},
    '30': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CE4R', 'Origin': 'E94'},
    '31': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOH', 'Origin': 'C94'},
    '32': {'alpha-i': '0.75', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'O2CM', 'Origin': 'C94'},
    '33': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HOS', 'Origin': 'C94'},
    '34': {'alpha-i': '1.00', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NR+', 'Origin': 'C94'},
    '35': {'alpha-i': '1.50', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OM', 'Origin': 'X94'},
    '36': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HNR+', 'Origin': 'C94'},
    '37': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CB', 'Origin': 'E94'},
    '38': {'alpha-i': '0.85', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NPYD', 'Origin': 'C94'},
    '39': {'alpha-i': '1.10', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NPYL', 'Origin': 'C94'},
    '40': {'alpha-i': '1.00', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NC=C', 'Origin': 'E94'},
    '41': {'alpha-i': '1.100', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CO2M', 'Origin': 'C94'},
    '42': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NSP', 'Origin': 'E94'},
    '43': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NSO2', 'Origin': 'E94'},
    '44': {'alpha-i': '3.00', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'STHI', 'Origin': 'C94'},
    '45': {'alpha-i': '1.150', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NO2', 'Origin': 'E94'},
    '46': {'alpha-i': '1.300', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'N=O', 'Origin': 'E94'},
    '47': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NAZT', 'Origin': 'X94'},
    '48': {'alpha-i': '1.200', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NSO', 'Origin': 'X94'},
    '49': {'alpha-i': '1.00', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'O+', 'Origin': 'X94'},
    '50': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HO+', 'Origin': 'C94'},
    '51': {'alpha-i': '0.400', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'O=+', 'Origin': 'E94'},
    '52': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HO=+', 'Origin': 'C94'},
    '53': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': '=N=', 'Origin': 'X94'},
    '54': {'alpha-i': '1.30', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'N+=C', 'Origin': 'C94'},
    '55': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NCN+', 'Origin': 'E94'},
    '56': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NGD+', 'Origin': 'E94'},
    '57': {'alpha-i': '1.000', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CNN+', 'Origin': 'E94'},
    '58': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NPD+', 'Origin': 'E94'},
    '59': {'alpha-i': '0.65', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OFUR', 'Origin': 'C94'},
    '60': {'alpha-i': '1.800', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'C%-', 'Origin': 'E94'},
    '61': {'alpha-i': '0.800', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NR%', 'Origin': 'E94'},
    '62': {'alpha-i': '1.300', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NM', 'Origin': 'X94'},
    '63': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C5A', 'Origin': 'E94'},
    '64': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C5B', 'Origin': 'E94'},
    '65': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5A', 'Origin': 'E94'},
    '66': {'alpha-i': '0.75', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5B', 'Origin': 'C94'},
    '67': {'alpha-i': '0.950', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N2OX', 'Origin': 'X94'},
    '68': {'alpha-i': '0.90', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N3OX', 'Origin': 'C94'},
    '69': {'alpha-i': '0.950', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'NPOX', 'Origin': 'C94'},
    '70': {'alpha-i': '0.87', 'N-i': '3.150', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'OH2', 'Origin': 'C94'},
    '71': {'alpha-i': '0.150', 'N-i': '0.800', 'A-i': '4.200', 'G-i': '1.209', 'DA': 'D', 'Symb': 'HS', 'Origin': 'C94'},
    '72': {'alpha-i': '4.000', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'SM', 'Origin': 'X94'},
    '73': {'alpha-i': '3.000', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': 'SMO2', 'Origin': 'X94'},
    '74': {'alpha-i': '3.000', 'N-i': '4.800', 'A-i': '3.320', 'G-i': '1.345', 'DA': '-', 'Symb': '=S=O', 'Origin': 'X94'},
    '75': {'alpha-i': '4.000', 'N-i': '4.500', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': '-P=C', 'Origin': 'X94'},
    '76': {'alpha-i': '1.200', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5M', 'Origin': 'X94'},
    '77': {'alpha-i': '1.500', 'N-i': '5.100', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'CLO4', 'Origin': 'X94'},
    '78': {'alpha-i': '1.350', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'C5', 'Origin': 'X94'},
    '79': {'alpha-i': '1.000', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5', 'Origin': 'X94'},
    '80': {'alpha-i': '1.000', 'N-i': '2.490', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'CIM+', 'Origin': 'C94'},
    '81': {'alpha-i': '0.80', 'N-i': '2.820', 'A-i': '3.890', 'G-i': '1.282', 'DA': '-', 'Symb': 'NIM+', 'Origin': 'C94'},
    '82': {'alpha-i': '0.950', 'N-i': '2.82', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'N5OX', 'Origin': 'X94'},
    '87': {'alpha-i': '0.45', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'FE+2', 'Origin': 'X94'},
    '88': {'alpha-i': '0.55', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'FE+3', 'Origin': 'X94'},
    '89': {'alpha-i': '1.4', 'N-i': '3.48', 'A-i': '3.890', 'G-i': '1.282', 'DA': 'A', 'Symb': 'F-', 'Origin': 'X94'},
    '90': {'alpha-i': '4.5', 'N-i': '5.100', 'A-i': '3.320', 'G-i': '1.345', 'DA': 'A', 'Symb': 'CL-', 'Origin': 'X94'},
    '91': {'alpha-i': '6.0', 'N-i': '6.000', 'A-i': '3.190', 'G-i': '1.359', 'DA': 'A', 'Symb': 'BR-', 'Origin': 'X94'},
    '92': {'alpha-i': '0.15', 'N-i': '2.', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'LI+', 'Origin': 'X94'},
    '93': {'alpha-i': '0.4', 'N-i': '3.5', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'NA+', 'Origin': 'X94'},
    '94': {'alpha-i': '1.0', 'N-i': '5.', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'K+', 'Origin': 'X94'},
    '95': {'alpha-i': '0.43', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'ZN+2', 'Origin': 'X94'},
    '96': {'alpha-i': '0.9', 'N-i': '5.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'CA+2', 'Origin': 'X94'},
    '97': {'alpha-i': '0.35', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'CU+1', 'Origin': 'X94'},
    '98': {'alpha-i': '0.40', 'N-i': '6.', 'A-i': '4.', 'G-i': '1.4', 'DA': '-', 'Symb': 'CU+2', 'Origin': 'X94'},
    '99': {'alpha-i': '0.35', 'N-i': '3.5', 'A-i': '4.', 'G-i': '1.3', 'DA': '-', 'Symb': 'MG+2', 'Origin': 'X94'},
}

# Amino acids and nucleic acid bases classification
molecule_class = {

    # Amino Acids
    "ALA": "hydrophobic",
    "ILE": "hydrophobic",
    "LEU": "hydrophobic",
    "VAL": "hydrophobic",
    "PHE": "hydrophobic",
    "PRO": "hydrophobic",
    "TRP": "hydrophobic",
    "MET": "hydrophobic",
    "GLY": "hydrophobic",
    "CYS": "polar charge: 0",
    "SER": "polar charge: 0",
    "THR": "polar charge: 0",
    "TYR": "polar charge: 0",
    "ASN": "polar charge: 0",
    "GLN": "polar charge: 0",
    "HIS": "polar charge: +",
    "LYS": "polar charge: +",
    "ARG": "polar charge: +",
    "ASP": "polar charge: -",
    "GLU": "polar charge: -",

    # Nucleic Acid Bases
    "A": "adenine base",
    "C": "cytosine base",
    "G": "guanine base",
    "T": "thymine base",   # DNA
    "U": "uracil base",    # RNA

    # Cofactors
    "MG" : "cofactor",
    "K"  : "cofactor",
    "CA" : "cofactor",
    "NA" : "cofactor",
    "HOH": "partial polar charge"
}

ionic_interactions = {

        "NA": ["F", "Cl", "Br", "I", "O"],
        "MG": ["F", "Cl", "Br", "I", "O"],
        "K" : ["F", "Cl", "Br", "I", "O"],
        "CA": ["F", "Cl", "Br", "I", "O"],
}

hydrogen_bond_acceptors = [
    "O",
    "N",
    "F"
    ]

# Hydrophobic residues and atoms for identifying hydrophobic interactions
hydrophobic_residues = [
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "TYR"
    ]

hydrophobic_atoms = [
    # Alanina
    "CB",
    # Valina
    "CB", "CG1", "CG2",
    # Isoleucina
    "CB", "CG1", "CG2", "CD1",
    # Leucina
    "CB", "CG", "CD1", "CD2",
    # Metionina
    "CB", "CG", "SD", "CE",
    # Prolina
    "CB", "CG", "CD",
    # Fenilalanina
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
    # Triptofano
    "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2",
    # Tirosina
    "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",

    # Átomos de carbono em grupos alquila e anéis alifáticos
    "C",
    # Nome comum para carbonos em anéis aromáticos
    "CA",
    # Outros nomes para carbonos em diferentes ambientes
    "CH",
    # Carbonos numerados, comuns em moléculas pequenas
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
    # Flúor
    "F",
    # Cloro
    "Cl",
    # Bromo
    "Br",
    # Iodo
    "I"
]

# Distance threshold for hydrophobic interactions
hydrophobic_distance_threshold = 4.0

 ######################################################################################################################


# Parser for PDB structure with consideration for protein structure, ligands, and cofactors
def parse_pdb_extended_v5(pdb_file):
    """
    Parses a PDB file to extract chains, cofactors, ligands, and atom details.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: A dictionary containing lists:
              - 'chains': List of ATOM details.
              - 'cofactors': List of HETATM details matching cofactor names.
              - 'ligands': List of HETATM details matching ligand names.
    """
    
    def extract_atom_details(line):
        """Helper function to extract atom details from a line."""
        return {
            'serial_number': int(line[6:11].strip()),
            'name': line[12:16].strip(),
            'alt_loc': line[16].strip(),
            'res_name': line[17:20].strip(),
            'chain_id': line[21].strip(),
            'res_seq': int(line[22:26].strip()),
            'icode': line[26].strip(),
            'coord': [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ],
            'occupancy': float(line[54:60].strip()),
            'temp_factor': float(line[60:66].strip()),
            'element': line[76:78].strip(),
            'charge': line[78:80].strip()
        }

    chains = []
    cofactors = []
    ligands = []

    cofactor_names = ["MG", "ZN", "CA", "K", "NA"]
    ligand_names = ["ACP", "TPS"]

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chains.append(extract_atom_details(line))
            elif line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                details = extract_atom_details(line)
                if residue_name in cofactor_names:
                    cofactors.append(details)
                elif residue_name in ligand_names:
                    ligands.append(details)

    return {
        'chains': chains,
        'cofactors': cofactors,
        'ligands': ligands
    }


def print_pdb_structure(pdb_dict):
    """
    Prints the PDB structure from the dictionary in an organized manner.

    Parameters:
        pdb_dict (dict): Dictionary containing the parsed PDB lists.

    Returns:
        None: Simply prints the PDB structured data.
    """

def format_line(atom_data, atom_type="ATOM  "):
    return (
        f"{atom_type:6s}{atom_data['serial_number']:5d} {atom_data['name']:<4s} {atom_data['alt_loc']:1s}{atom_data['res_name']:<3s} "
        f"{atom_data['chain_id']:1s}{atom_data['res_seq']:4d}{atom_data['icode']:1s}   "
        f"{atom_data['coord'][0]:8.3f}{atom_data['coord'][1]:8.3f}{atom_data['coord'][2]:8.3f}"
        f"{atom_data['occupancy']:6.2f}{atom_data['temp_factor']:6.2f}          "
        f"{atom_data['element']:^2s}{atom_data['charge']:2s}\n"
    )

def print_pdb_structure(pdb_dict):
    # Print the chains first
    for atom in pdb_dict['chains']:
        print(format_line(atom), end='')

    # Print cofactors
    print("TER")
    for atom in pdb_dict['cofactors']:
        print(format_line(atom, "HETATM"), end='')

    # Print ligands
    print("TER")
    for atom in pdb_dict['ligands']:
        print(format_line(atom, "HETATM"), end='')

    print("END")


def map_to_molecular_group(atom_name, residue_name):

    residue_name = residue_name.upper()  # Garante que o nome do resíduo esteja em maiúsculas
    atom_name = atom_name.upper()  # Similarmente, para o nome do átomo

    # Lista de aminoácidos com carbonos laterais específicos
    lateral_carbon_aminoacids = ["ALA", "VAL", "LEU", "ILE", "PRO", "MET", "CYS", 
                                 "PHE", "TYR", "THR", "LYS", "ARG", "GLU", "ASP"]
    
    if atom_name.startswith("C"):

        # Verificar condição do grupo alquila
        if residue_name == "VAL" and atom_name in ["CG1", "CG2"]:
            return "CR"  
        elif residue_name == "LEU" and atom_name in ["CG", "CD1", "CD2"]:
            return "CR"
        elif residue_name == "ILE" and atom_name in ["CG1", "CG2", "CD1"]:
            return "CR"
        elif residue_name == "MET" and atom_name in ["CG", "SD", "CE"]:
            return "CR"  
        # Verificando a condição de carbono alfa
        elif atom_name == "CA" and residue_name != "PRO":
            return "HC"
        # Verificando condição de carbonos laterais
        elif atom_name.startswith("C") and residue_name in lateral_carbon_aminoacids:
            return "HC"
        # Condição para os carbonos aromáticos
        if residue_name == "PHE" and atom_name.startswith("CZ"):  # Ajuste a nomeação do átomo conforme necessário
            return "CB"
        elif residue_name == "TYR" and atom_name.startswith("CZ"):  # Ajuste a nomeação do átomo conforme necessário
            return "CB"
        elif residue_name == "TRP" and atom_name in ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]:  # A nomeação dos átomos aqui é apenas uma suposição, ajuste conforme necessário
            return "CB"
        elif residue_name in ["ASP", "GLU"] and atom_name in ["CG", "CD"]:  # CG e CD são átomos de carbono em Asp e Glu, respectivamente, que fazem parte do grupo carboxila
            return "CO2M"
        # Condição para o carbono no grupo guanidínio da Arginina
        elif residue_name == "ARG" and atom_name == "CG":
            return "CGD"
        else: 
            return "CR"

    if atom_name.startswith("O"):

        # Condição de "O=C"
        if atom_name in ["OD1", "OE1", "OE2", "O"]:
            return "O=C"
        # Termino carboxila de todos os aminoácidos
        elif atom_name in ["OD2", "OG", "OG1", "OH"]:
            return "OR"
        # Condição de Serina
        elif residue_name == "SER" and atom_name == "OG":
            return "OR"
        # Condição de Treonina
        elif residue_name == "THR" and atom_name == "OG1":
            return "OR"
        # Condição de Tirosina
        elif residue_name == "TYR" and atom_name == "OH":
            return "OR"
        # Condição para Ácido aspártico e Ácido glutâmico
        elif residue_name in ["ASP", "GLU"] and atom_name in ["OD1", "OD2", "OE1", "OE2"]:
            return "OR"
        # Condição para Asparagina e Glutamina
        elif residue_name in ["ASN", "GLN"] and atom_name in ["OD1", "OE1"]:
            return "OR"
        # Condição para Ácido aspártico e Ácido glutâmico
        elif residue_name in ["ASP", "GLU"] and atom_name in ["OD1", "OD2", "OE1", "OE2"]:
            return "O2CM"
        else: 
            return "OR"

    if atom_name.startswith("N"):

        # Verifica o nitrogênio do grupo amina terminal para todos os aminoácidos
        if atom_name == "N":
            return "NR"
        # Condição específica para Lisina
        elif residue_name == "LYS" and atom_name in ["NZ"]:
            return "NR"
        # Condição específica para Arginina
        elif residue_name == "ARG" and atom_name in ["NE", "NH1", "NH2"]:
            return "NGD+"
        # Condição específica para Histidina
        elif residue_name == "HIS" and atom_name in ["ND1", "NE2"]:
            return "NPYD"
        # Condição para grupo Amido (N-C=O)
        elif (residue_name == "ASN" and atom_name == "ND2") or \
        (residue_name == "GLN" and atom_name == "NE2"):
            return "NC=O"
        # Condição para o nitrogênio no anel pirrol do triptofano
        elif residue_name == "TRP" and atom_name == "NE1":  # NE1 é o nitrogênio no anel de pirrol do triptofano
            return "NPYL"
        elif residue_name == "PRO" and atom_name == "N":  # N é o nitrogênio da prolina que forma a ligação imina
            return "'NC=C"
        # Condição para o nitrogênio com ligação tripla em uma variante modificada da Cisteína
        elif residue_name == "CYS" and atom_name == "NSP":  # NSP seria um nome hipotético para esse nitrogênio com ligação tripla
            return "NSP"
        else: 
            return "NR"

    if atom_name.startswith("H"):

        # Checa se o átomo é um hidrogênio ligado a um nitrogênio amida em um peptídeo/proteína (Excluindo resíduos terminais). 
        if atom_name in ["HN"] or atom_name.startswith("H") and atom_name[1:].isdigit():
            return "HNCO"
        elif residue_name in ["SER", "THR"] and "HO" in atom_name: 
            return "HOR"
        elif residue_name == "LYS" and "HNZ" in atom_name:
            return "HNR"
        elif residue_name in ["ASP", "GLU"] and "HOCO" in atom_name:
            return "HOCO"
        elif residue_name in ["SER", "THR"] and "HO" in atom_name: 
            return "HOH"
        else:
            return "HR"
    
    if atom_name.startswith("S"):

         # Condição para grupos tiol e tioéter
        if (residue_name == "MET" and atom_name == "SD"):
            return "S"
        # Condição para o átomo de enxofre no grupo tiol da Cisteína
        elif residue_name == "CYS" and atom_name == "SG":
            return "HS"
             

def find_molecule(pdb_dict, molecule_name):
    """
    Searches for a molecule in the parsed PDB dictionary based on its name.

    Parameters:
        pdb_dict (dict): Dictionary containing parsed PDB data.
        molecule_name (str): Name of the molecule to search for.

    Returns:
        tuple: (Name of the molecule, residue sequence number) or None if not found.
    """

    chain_select = None

    # Check if the chain_select was provided as an argument
    if len(sys.argv) > 4:
        chain_select = str(sys.argv[4])

    # Combine all atom lists
    all_atoms = pdb_dict['chains'] + pdb_dict['cofactors'] + pdb_dict['ligands']

    # Check if any atom with the desired molecule name is present
    for atom in all_atoms:
        if atom['res_name'] == molecule_name:
            if chain_select is None or atom['chain_id'] == chain_select:
                return (atom['res_name'], atom['res_seq'], atom['chain_id'])
            else:
                atom['chain_id'] = chain_select
                return (atom['res_name'], atom['res_seq'], atom['chain_id'])



def verify_near_residues(input_pdb, ligand_residue, treshold_distance):

    # Extract ligand atoms based on ligand_residue tuple
    ligand_atoms = [atom for atom in input_pdb['ligands'] 
                    if atom['res_name'] == ligand_residue[0] 
                    and atom['res_seq'] == ligand_residue[1] 
                    and atom['chain_id'] == ligand_residue[2]]

    # List to store the next residues and the interaction
    near_residues_info = []

    # Check all atoms of all residues
    all_atoms = input_pdb['chains'] + input_pdb['cofactors'] + input_pdb['ligands']

    for atom in all_atoms:
        if (atom['res_name'] != ligand_residue[0] or 
           atom['res_seq'] != ligand_residue[1] or 
           atom['chain_id'] != ligand_residue[2]):  # Ignore specified ligand
            min_distance = float('inf')
            interacting_atoms = None

            for ligand_atom in ligand_atoms:
                # Calculate the euclidean distance manually
                distance = ((atom['coord'][0] - ligand_atom['coord'][0]) ** 2 +
                            (atom['coord'][1] - ligand_atom['coord'][1]) ** 2 +
                            (atom['coord'][2] - ligand_atom['coord'][2]) ** 2) ** 0.5

                if distance < min_distance:
                    min_distance = distance
                    interacting_atoms = (atom, ligand_atom)

            if min_distance <= treshold_distance:
                info = {
                    'molecule_name': atom['res_name'],
                    'molecule_number': atom['res_seq'],
                    'chain': atom['chain_id'],
                    'distance': min_distance,
                    'molecule_atom': atom['name'],
                    'ligand_atom': interacting_atoms[1]['name']
                }
                near_residues_info.append(info)

    # Sort the near_residues_info list based on the distance
    near_residues_info.sort(key=lambda x: x['distance'])

    return near_residues_info



def set_output(output_name, near_residues_dict, ligand_residue_tuple):
    ligand_name, ligand_num, ligand_chain = ligand_residue_tuple

    interacting_molecules_count = 0  # contador para moléculas interagindo

    with open(output_name, 'w', newline='') as file:
        writer = csv.writer(file)
        columns = ["Molecule", "Classification", "Number", "Chain", "Nearby atoms", "Distance(Å)", "Interaction"]
        writer.writerow(columns)

        print("{:^20} {:^30} {:^10} {:^5} {:^20} {:^10} {:^20}".format(*columns))

        for entry in near_residues_dict:
            aa_name = entry['molecule_name']
            aa_num = entry['molecule_number']
            chain_id = entry['chain']
            distance = entry['distance']
            molecule_atom = entry['molecule_atom']
            ligand_atom = entry['ligand_atom']

            # Use molecule_class para obter a classificação, se disponível, caso contrário, defina como 'unknown'
            aa_class = molecule_class.get(aa_name, "unknown")

            atom1_str = molecule_atom + "(" + aa_name + ")"
            atom2_str = ligand_atom + "(" + ligand_name + ")"
            nearby_atoms_str = "{:<10}-{:>10}".format(atom1_str, atom2_str)

            # Supondo que sua função is_interaction funcione assim
            probable_interaction = is_interaction(molecule_atom, ligand_atom, aa_name, distance)

            # Incrementa o contador se houver interação
            if probable_interaction:
                interacting_molecules_count += 1

            writer.writerow([aa_name, aa_class, aa_num, chain_id, nearby_atoms_str, round(distance, 2), probable_interaction])

            print("{:^20} {:^30} {:^10} {:^5} {:^20} {:^10.2f} {:^20}".format(
                aa_name, aa_class, aa_num, chain_id, nearby_atoms_str, distance, probable_interaction))

    print("\nTotal number of interacting molecules:", interacting_molecules_count)
    print("\n")
    return f"Successfully processed and saved! Total interactions: {interacting_molecules_count}"



# Define a function to check for potential hydrogen bonds
def is_interaction(atom1_name, atom2_name, residue_name, distance):

    # Check for hydrophobic interactions
    if residue_name in hydrophobic_residues:
        if atom1_name in hydrophobic_atoms or atom2_name in hydrophobic_atoms:
            if distance <= hydrophobic_distance_threshold:
                return "Hydrophobic"
            
    # Check for ionic interaction
    if distance < 4:
        if residue_name in ionic_interactions:
            if atom1_name.startswith(tuple(ionic_interactions.get(residue_name, []))) or \
                atom2_name.startswith(tuple(ionic_interactions.get(residue_name, []))):
                return "Ionic"
            
    if distance < 3:
        # Check for hydrogen bond
        if atom1_name.startswith(tuple(hydrogen_bond_acceptors)) and \
           atom2_name.startswith(tuple(hydrogen_bond_acceptors)):
            return "Hydrogen bond"

    # Assuming lennard_jones_potential is defined elsewhere in your code
    v_lj = lennard_jones_potential(atom1_name, atom2_name, residue_name, distance)

    # Check for potential van der Waals interaction based on Lennard-Jones potential
    if v_lj < -0.1:
        return "van der Waals"

    return "Non-specific"


def lennard_jones_potential(atom1_name, atom2_name, residue, r):
    """Calculate Lennard-Jones potential between two atoms based on MMFF types."""
    # Get MMFF types for the atoms
    mapped_group1 = map_to_molecular_group(atom1_name, residue)
    mapped_group2 = map_to_molecular_group(atom2_name, residue)
    if not mapped_group1 or not mapped_group2: 
        return 0  # or handle this case as required

    type1 = molecular_group_dict[mapped_group1]['PRIMARY MMF TYPE']
    type2 = molecular_group_dict[mapped_group2]['PRIMARY MMF TYPE']

    # Get epsilon and sigma values for the atoms
    epsilon1 = float(vdw_dict[type1]['alpha-i'])
    epsilon2 = float(vdw_dict[type2]['alpha-i'])

    sigma1 = float(vdw_dict[type1]['N-i'])
    sigma2 = float(vdw_dict[type2]['N-i'])

    # Combine the epsilon and sigma values
    epsilon_combined = (epsilon1 * epsilon2) ** 0.5
    sigma_combined = (sigma1 + sigma2) / 2.0

    # Calculate the Lennard-Jones potential
    return 4 * epsilon_combined * ((sigma_combined / r)**12 - (sigma_combined / r)**6)



def format_output(ligand_residue_tuple, near_residues_dict):
    ligand, ligand_number = ligand_residue_tuple
    formatted_results = []

    for entry in near_residues_dict:
        line = (f"LIGAND RESIDUE: {ligand} {ligand_number} | "
                f"NEAR RESIDUE: {entry['molecule_name']} {entry['molecule_number']} {entry['chain']} "
                f"ATOM: {entry['molecule_atom']} OF {entry['molecule_name']} "
                f"DISTANCE TO {entry['ligand_atom']} OF {ligand} = {entry['distance']:.2f} A")
        
        formatted_results.append(line)

    return formatted_results


if __name__ == "__main__":

    # Arquivos de entrada e saida a serem fornecidos
    input_pdb      = sys.argv[1]    # example.pdb
    input_molecule = sys.argv[2]    # ATP
    output_name    = sys.argv[3]    # ATP_OUT (csv)

    # Distance from the selected molecule
    treshold_distance = 4.0
    input_pdb = parse_pdb_extended_v5(input_pdb)

    # executa e salva o resultados para a classificacao dos contatos
    ligand_residue = find_molecule(input_pdb, input_molecule)
    near_residues  = verify_near_residues(input_pdb, ligand_residue, treshold_distance)
    set_output(output_name, near_residues, ligand_residue)


# verificar artigos sobre cinases + IA