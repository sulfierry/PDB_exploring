# Check possible contacts between structure and ligand
# python script.py ref.pdb LIG PDB_LIG.csv CHAIN(optional)
# GitHub: github.com/sulfierry/

import csv
import sys
import math


# Variaveis globais ################################################################################################

aminoacid_group_dict = {

     # Carbonos da cadeia lateral de valina, leucina, isoleucina e metionina são todos carbonos alquílicos
     # Valina (VAL): CG1, CG2
     # Leucina (LEU): CG, CD1, CD2
     # Isoleucina (ILE): CG1, CG2, CD1
     # Metionina (MET): CG, SD, CE
    'CR': {'PRIMARY MMF TYPE': '1', 'DEFAULT TYPES': ['1', '1', '1', '0'], 'DEFINITION': 'ALKYL CARBON'},
    

    # Alanina (Ala, A) - Grupo lateral: CH3 . Valina (Val, V) - Grupo lateral: CH(CH3)2. Leucina (Leu, L) - Grupo lateral: CH2CH(CH3)2
    # Isoleucina (Ile, I) - Grupo lateral: CH(CH3)CH2CH3. Proleína (Pro, P) - Grupo lateral: Anel alifático que se conecta ao nitrogênio e ao carbono alfa.
    # Metionina (Met, M) - Grupo lateral: CH2CH2SCH3. Cisteína (Cys, C) - Grupo lateral: CH2SH Fenilalanina (Phe, F) - Grupo lateral: CH2 com anel aromático.
    # Tirosina (Tyr, Y) - Grupo lateral: CH2 com anel aromático e OH. Treonina (Thr, T) - Grupo lateral: CH(OH)CH3. Lisina (Lys, K) - Grupo lateral: CH2 (4 vezes) terminando em NH3+.
    # Arginina (Arg, R) - Grupo lateral: CH2 (3 vezes) seguido de um grupo guanidino. Glutamato (Glu, E) e Aspartato (Asp, D) - Estes possuem uma parte de sua cadeia lateral que é CH2.
    # Observe que o esqueleto comum de aminoácidos (parte fora do grupo lateral) também possui ligações H-C.  
    # O carbono alfa, que está ligado diretamente ao grupo amino e ao grupo carboxila, está ligado a um hidrogênio em todos os aminoácidos exceto a prolina.
    'HC': {'PRIMARY MMF TYPE': '5', 'DEFAULT TYPES': ['5', '5', '5', '0'], 'DEFINITION': 'H-C'},

    # O termo "O-CSP3" se refere a um átomo de oxigênio ligado a um átomo de carbono sp3 (carbono tetraédrico/saturado). 
    # Nos aminoácidos padrão, o "O-CSP3" é comumente encontrado no grupo carboxila (COOH) e, para alguns aminoácidos, em grupos laterais específicos.
    # Todos os aminoácidos: Na terminação carboxila (COOH). Portanto, todos os aminoácidos padrão terão pelo menos um O-CSP3 no carbono da carboxila. \
    # Serina (Ser, S): Grupo lateral: CH2OH. O oxigênio do grupo hidroxila (OH) está ligado a um carbono sp3. 
    # Treonina (Thr, T): Grupo lateral: CH(OH)CH3. O oxigênio do grupo hidroxila (OH) está ligado a um carbono sp3. 
    # Tirosina (Tyr, Y): Grupo lateral: CH2 com anel aromático e OH. O oxigênio do grupo hidroxila (OH) está ligado a um carbono sp3. 
    # Ácido aspártico (Asp, D) e Ácido glutâmico (Glu, E): Ambos têm um segundo O-CSP3 no grupo carboxilato de suas cadeias laterais. 
    # Asparagina (Asn, N) e Glutamina (Gln, Q): Ambos têm um O-CSP3 em suas cadeias laterais ligado ao carbono do grupo amida.
    'OR': {'PRIMARY MMF TYPE': '6', 'DEFAULT TYPES': ['6', '6', '6', '0'], 'DEFINITION': 'O-CSP3'},

    # O termo "O=C" refere-se ao oxigênio de uma ligação dupla carbonílica, que é encontrado no grupo carboxila (COOH) e em outros grupos funcionais de aminoácidos.
    # O grupo carboxila (COOH) presente em todos os aminoácidos padrão tem um oxigênio carbonílico. Portanto, todos os aminoácidos padrão terão pelo menos um "O=C" devido à sua terminação carboxila.
    # Ácido aspártico (Asp, D) e Ácido glutâmico (Glu, E): Além do grupo carboxila, eles também têm um grupo carboxilato adicional em suas cadeias laterais, fornecendo um segundo "O=C".
    # Asparagina (Asn, N) e Glutamina (Gln, Q): Ambos têm um grupo amida em suas cadeias laterais, que contém um "O=C". Tirosina (Tyr, Y): Contém um grupo fenol em sua cadeia lateral, 
    # e embora o anel benzênico contenha oxigênios, eles não estão na forma "O=C". No entanto, estou mencionando a tirosina para esclarecer que, embora contenha oxigênio, não tem o grupo "O=C" na cadeia lateral.
    # Cisteína (Cys, C): Possui um grupo tiol, mas, como se refere a enxofre, não é "O=C". Mais uma vez, estou mencionando a cisteína apenas para esclarecimento.
    'O=C': {'PRIMARY MMF TYPE': '7', 'DEFAULT TYPES': ['7', '7', '6', '0'], 'DEFINITION': 'O=C, GENERIC'},
    
    # O termo "AMINE N" refere-se ao nitrogênio presente em aminas. Nos aminoácidos, o nitrogênio da amina é tipicamente encontrado no grupo amina (NH2, NH ou N).
    # Todos os aminoácidos: Todos os 20 aminoácidos padrão têm um grupo amina na terminação N-terminal. Portanto, eles todos têm pelo menos um nitrogênio da amina.
    # Lisina (Lys, K): Além do grupo amina no N-terminal, a lisina possui um segundo grupo amina em sua cadeia lateral. Arginina (Arg, R): A arginina tem um grupo 
    # guanidina em sua cadeia lateral, que contém múltiplos nitrogênios. No entanto, nem todos esses nitrogênios são equivalentes a nitrogênios de amina típicos, 
    # mas estão presentes na estrutura geral do grupo guanidina. Histidina (His, H): A histidina tem um anel imidazol em sua cadeia lateral que contém dois nitrogênios. 
    # Novamente, embora estes não sejam nitrogênios de aminas típicas, eles são nitrogênios em um contexto de heterociclo.
    'NR': {'PRIMARY MMF TYPE': '8', 'DEFAULT TYPES': ['8', '8', '8', '0'], 'DEFINITION': 'AMINE N'},


    # Os amidos (N-C=O) são encontrados nos aminoácidos que possuem cadeias laterais contendo o grupo funcional conamida.
    # Asparagina (Asn, N): Cadeia lateral: CH₂CONH₂. Posição atômica do amido: A nitrogênio da amida é denominado ND2, e o carbono da amida é denominado CG. 
    # Glutamina (Gln, Q): Cadeia lateral: CH₂CH₂CONH₂. Posição atômica do amido: O nitrogênio da amida é denominado NE2, e o carbono da amida é denominado CD. 
    # Além disso, a ligação peptídica que liga aminoácidos em uma proteína é, por definição, uma ligação amida, mas ela é formada entre o grupo carboxílico de 
    # um aminoácido e o grupo amina do próximo, e não é parte da cadeia lateral de qualquer aminoácido.
    'NC=O': {'PRIMARY MMF TYPE': '10', 'DEFAULT TYPES': ['10', '10', '8', '0'], 'DEFINITION': 'N-C=O, AMIDES'},

    # Os aminoácidos que contêm grupos tiol (–SH) e sulfeto (ou dissulfeto, –S–S–) são:
    # Cisteína (Cys, C): Contém um grupo tiol (–SH). Quando duas moléculas de cisteína estão próximas uma da outra em uma proteína, 
    # elas podem formar uma ligação dissulfeto (–S–S–), que é covalente e ajuda a estabilizar a estrutura tridimensional da proteína. 
    # Metionina (Met, M): Não possui um grupo tiol, mas possui um átomo de enxofre em sua cadeia lateral na forma de um grupo tioéter. 
    # Embora não forme ligações dissulfeto como a cisteína, a metionina é relevante quando se fala de aminoácidos contendo enxofre.
    'S': {'PRIMARY MMF TYPE': '15', 'DEFAULT TYPES': ['15', '15', '15', '0'], 'DEFINITION': 'THIOL, SULFIDE'},


    # O aminoácido padrão metionina não possui um enxofre na forma de sulfona em sua estrutura natural. 
    # No entanto, quando a metionina é exposta a agentes oxidantes mais fortes ou por um período mais longo, pode ocorrer uma oxidação adicional 
    # do sulfoxido de metionina para formar sulfona de metionina. Em outras palavras, a metionina pode ser oxidada a metionina sulfoxido e, em 
    # seguida, a metionina sulfoxido pode ser oxidada adicionalmente para formar sulfona de metionina.
    'SO2': {'PRIMARY MMF TYPE': '18', 'DEFAULT TYPES': ['18', '18', '15', '0'], 'DEFINITION': 'SULFONE S'},

    # Em relação aos aminoácidos padrão encontrados nas proteínas, os que contêm um grupo hidroxila (H-O) e, portanto, podem ser considerados alcoóis são: 
    # Serina (Ser, S): Este aminoácido possui um grupo lateral -CH2OH. Treonina (Thr, T): Contém um grupo lateral -CHOHCH3.
    'HOR': {'PRIMARY MMF TYPE': '21', 'DEFAULT TYPES': ['21', '21', '5', '0'], 'DEFINITION': 'H-O, ALCOHOLS'},

    # Aminas primárias são caracterizadas pela presença de um grupo amino (-NH2). Os aminoácidos que possuem um grupo amino (em adição ao grupo amino terminal, que todos os aminoácidos têm) são: 
    # Lisina (Lys, K): Possui uma cadeia lateral que termina em um grupo -NH2. A posição atômica em uma molécula típica de lisina para este nitrogênio é Nζ (ou, em alguns formatos de notação, NZ). 
    # Arginina (Arg, R): Sua cadeia lateral contém três grupos nitrogênio, mas eles são parte de um grupo guanidino. Assim, embora a arginina contenha nitrogênios, eles não estão na forma de aminas primárias típicas.
    'HNR': {'PRIMARY MMF TYPE': '23', 'DEFAULT TYPES': ['23', '23', '5', '0'], 'DEFINITION': 'H-N, AMINES'},

    # Quando se refere a "H-O, ACIDS", parece estar se referindo aos grupos hidroxila que fazem parte de ácidos carboxílicos. No contexto dos aminoácidos, 
    # todos eles possuem um grupo carboxílico (-COOH) na extremidade C-terminal, mas este grupo é geralmente desprotonado em pH fisiológico, tornando-se -COO^-. 
    # Dentre os aminoácidos padrão, temos: Ácido aspártico (Asp, D): Possui um grupo carboxílico adicional em sua cadeia lateral. Em pH fisiológico, este grupo 
    # geralmente está na forma desprotonada, -COO^-. Ácido glutâmico (Glu, E): Similar ao ácido aspártico, também tem um grupo carboxílico extra em sua cadeia 
    # lateral que, em pH fisiológico, estará na forma -COO^-.
    'HOCO': {'PRIMARY MMF TYPE': '24', 'DEFAULT TYPES': ['24', '24', '5', '0'], 'DEFINITION': 'H-O, ACIDS'},


    # Todos os aminoácidos em peptídeos e proteínas, exceto os resíduos terminais, apresentarão um átomo de nitrogênio ligado por uma ligação 
    # simples a um átomo de hidrogênio (H-N) em sua estrutura amida.
    'HNCO': {'PRIMARY MMF TYPE': '28', 'DEFAULT TYPES': ['28', '28', '5', '0'], 'DEFINITION': 'H-N, AMIDES'},

    #  A serina e a treonina, dois dos 20 aminoácidos padrão, têm um grupo hidroxila (-OH) como parte de suas cadeias laterais. Serina (Ser, S) tem a cadeia lateral: -CH₂-OH 
    # Treonina (Thr, T) tem a cadeia lateral: -CH(OH)-CH₃
    'HOH': {'PRIMARY MMF TYPE': '31', 'DEFAULT TYPES': ['31', '31', '31', '0'], 'DEFINITION': 'H-OH'},

    # O grupo "O, CARBOXYLATE ANION" refere-se ao oxigênio encontrado no ânion carboxilato (COO⁻).
    # Dos 20 aminoácidos padrão, somente dois são aminoácidos ácidos (também chamados de aminoácidos dicarboxílicos) porque possuem um grupo carboxilato adicional 
    # em suas cadeias laterais quando estão no estado fisiológico (pH ~7.4). Esses aminoácidos são: 
    # Ácido Aspártico (Asp, D): Sua cadeia lateral é -CH₂-COOH, mas em pH fisiológico, a forma predominante é -CH₂-COO⁻, onde o grupo carboxilato está desprotonado. 
    # Ácido Glutâmico (Glu, E): Sua cadeia lateral é -CH₂-CH₂-COOH, e em pH fisiológico, a forma predominante é -CH₂-CH₂-COO⁻, com o grupo carboxilato desprotonado. 
    # Assim, em ambientes de pH ao redor do neutro, esses aminoácidos terão um grupo carboxilato anionico (COO⁻) em suas cadeias laterais.
    'O2CM': {'PRIMARY MMF TYPE': '32', 'DEFAULT TYPES': ['32', '7', '6', '0'], 'DEFINITION': 'O, CARBOXYLATE ANION'},

    # Dentre os 20 aminoácidos padrão, a Lisina (Lys, K) é o aminoácido que possui potencial para ter um nitrogênio quaternário em sua cadeia lateral, 
    # devido à sua terminação em amina primária. Em um ambiente fisiológico, com pH ao redor de 7,4, a terminação amina da lisina é frequentemente protonada, 
    # adotando uma carga positiva, mas não é um nitrogênio quaternário verdadeiro porque só possui três substituintes.
    'NR+': {'PRIMARY MMF TYPE': '34', 'DEFAULT TYPES': ['34', '8', '8', '0'], 'DEFINITION': 'N+, QUATERNARY N'},

    # O grupo funcional "H-N+" corresponde a um átomo de hidrogênio ligado a um átomo de nitrogênio positivamente carregado, formando um cátion. 
    # Esse tipo de grupo funcional é encontrado em aminoácidos que possuem um grupo amina protonado. Os aminoácidos que podem apresentar o grupo H-N+ 
    # são aqueles cujos grupos amina passaram por uma protonação para formar íons aminium.
    # Lisina (Lys): A lisina possui uma cadeia lateral que contém um grupo amina (NH2) na extremidade da cadeia. Quando a lisina ganha um próton (H+) nesse grupo amina, 
    # ela forma o íon aminium (H3N+), que possui o grupo funcional H-N+.  Arginina (Arg): A arginina também possui uma cadeia lateral que contém um grupo guanidino (NH=C(NH2)2) 
    # que é básico e pode facilmente ganhar um próton para formar o íon aminium (H3N+), criando o grupo funcional H-N+.
    'HNR+': {'PRIMARY MMF TYPE': '36', 'DEFAULT TYPES': ['36', '36', '5', '0'], 'DEFINITION': 'H-N+'}, # Deixei o valor NR

    # Os aminoácidos que apresentam um grupo funcional "AROMATIC C" são aqueles que possuem um anel aromático de carbono em sua estrutura
    # Fenilalanina (Phe): A fenilalanina tem uma cadeia lateral que consiste em um anel benzênico, que é um exemplo clássico de anel aromático. 
    # Tirosina (Tyr): A tirosina também possui um anel benzênico em sua cadeia lateral, semelhante ao da fenilalanina. Além disso, a tirosina contém um grupo hidroxila (-OH) na posição para.
    # Triptofano (Trp): O triptofano tem uma cadeia lateral que contém um anel indol, que é um tipo de anel aromático derivado do benzênico. Esse anel contém um átomo de nitrogênio e é característico do triptofano.
    'CB': {'PRIMARY MMF TYPE': '37', 'DEFAULT TYPES': ['37', '2', '1', '0'], 'DEFINITION': 'AROMATIC C'},
    
    # A piridina é um anel heterocíclico de seis membros, contendo um átomo de nitrogênio em uma das posições do anel.
    # O aminoácido que contém um anel de piridina em sua estrutura é a histidina (His). 
    'NPYD': {'PRIMARY MMF TYPE': '38', 'DEFAULT TYPES': ['38', '9', '8', '0'], 'DEFINITION': 'AROMATIC N, PYRIDINE'},

    # O pirrol é um anel heterocíclico de cinco membros, contendo um átomo de nitrogênio em uma das posições do anel.
    # O aminoácido que contém um anel de pirrol em sua estrutura é o triptofano (Trp). O triptofano é conhecido por 
    # possuir um anel indol em sua cadeia lateral, que é um tipo de anel aromático derivado do benzênico. 
    # O indol é um derivado do pirrol, contendo um grupo nitrogênio em uma das posições do anel, assim caracterizando o "AROMATIC N, PYRROLE".
    'NPYL': {'PRIMARY MMF TYPE': '39', 'DEFAULT TYPES': ['39', '10', '8', '0'], 'DEFINITION': 'AROMATIC N, PYRROLE'},

    # O grupo funcional "N-C=C (DELOC LP)" se refere a um átomo de nitrogênio (N) ligado a um átomo de carbono (C) com uma ligação dupla, 
    # onde a carga da ligação dupla é delocalizada, ou seja, o par de elétrons da ligação está compartilhado entre os dois átomos. 
    # Esse tipo de ligação é chamado de ligação dupla imina (ou ligação dupla C=N). O aminoácido que contém o grupo funcional "N-C=C (DELOC LP)" é a prolina (Pro).
    # A prolina tem uma estrutura única na qual a sua cadeia lateral forma uma ligação imina com o grupo amino do esqueleto do aminoácido, resultando em uma conformação cíclica.
    'NC=C': {'PRIMARY MMF TYPE': '40', 'DEFAULT TYPES': ['40', '10', '8', '0'], 'DEFINITION': 'N-C=C (DELOC LP)'},

    # O grupo funcional "C IN CO2- ANION" refere-se a um átomo de carbono (C) presente no íon carboxilato (CO2-). 
    # O íon carboxilato é formado quando um grupo carboxila (COOH) perde um próton (H+) e se torna carregado negativamente.
    # Os aminoácidos com grupos carboxila são: ácido aspártico (Asp) e ácido glutâmico (Glu).
    # Estes com cadeias laterais carregadas ou polares incluem o grupo carboxila, que pode perder um próton para formar um íon carboxilato.
    'CO2M': {'PRIMARY MMF TYPE': '41', 'DEFAULT TYPES': ['41', '3', '1', '0'], 'DEFINITION': 'C IN CO2- ANION'},

    # O grupo funcional "N TRIPLE BONDED" se refere a um átomo de nitrogênio (N) ligado a outro átomo por uma ligação tripla (triple bond). Essa ligação tripla pode ser chamada de ligação nitrila (C≡N).
    # O único aminoácido que contém uma ligação tripla nitrila é a cisteína (Cys). A cisteína normalmente forma pontes dissulfeto com outras cisteínas, mas sob certas condições, pode ser modificada para formar um grupo nitrila na cadeia lateral.
    'NSP': {'PRIMARY MMF TYPE': '42', 'DEFAULT TYPES': ['42', '42', '8', '0'], 'DEFINITION': 'N TRIPLE BONDED'},

    # O termo "GUANIDINIUM N; Q=1/3" refere-se a um átomo de nitrogênio (N) em um grupo funcional guanidínio, onde o íon guanidínio possui uma carga positiva (Q=1/3). 
    # O grupo guanidínio consiste em três átomos de nitrogênio ligados a um átomo de carbono, e é encontrado em aminoácidos como a arginina. 
    # Na arginina, um dos aminoácidos padrão, a cadeia lateral contém um grupo guanidínio, que é importante para interações em proteínas e em várias reações bioquímicas. 
    'NGD+': {'PRIMARY MMF TYPE': '56', 'DEFAULT TYPES': ['56', '10', '8', '0'], 'DEFINITION': 'GUANIDINIUM N; Q=1/3'},

    # O termo "GUANIDINIUM CARBON" se refere a um átomo de carbono (C) em um grupo funcional guanidínio. O grupo guanidínio é composto por três átomos de nitrogênio (N) 
    # e um átomo de carbono (C) ligados entre si de forma específica. A arginina é um aminoácido que contém o grupo funcional guanidínio, onde um dos átomos de carbono 
    # está envolvido na estrutura do grupo.
    'CGD+': {'PRIMARY MMF TYPE': '57', 'DEFAULT TYPES': ['57', '2', '1', '0'], 'DEFINITION': 'GUANIDINIUM CARBON'},

    # A cisteína é um aminoácido que contém um grupo tiol (-SH), que pode formar ligações de dissulfeto com outros grupos tiol em outras moléculas de cisteína. 
    'HS': {'PRIMARY MMF TYPE': '71', 'DEFAULT TYPES': ['71', '5', '5', '0'], 'DEFINITION': 'H-S'},
    
    'MG+2': {'PRIMARY MMF TYPE': '99', 'DEFAULT TYPES': ['99', '99', '99', '99'], 'DEFINITION': 'DIPOSITIVE MAGNESIUM CATION'},

}

aminoacid_vdw_dict = {
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
    "HOH": "partial polar charge",

    # ligands
    "TPS": "ligand",
    "ACP": "ligand",
    "TMP": "ligand",
    "TPP": "ligand",
    "ANP": "ligand",   
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
    "F",
    "H"
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
    "H",
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

ionic_atoms = {
    'positive': ['NH1', 'NH2', 'NZ', 'ND1', 'NE2'],  # Átomos carregados positivamente
    'negative': ['OD1', 'OD2', 'OE1', 'OE2']  # Átomos carregados negativamente
}


# Distance threshold for hydrophobic interactions
hydrophobic_distance_threshold = 4.0

 ######################################################################################################################


# Parser for PDB structure with consideration for protein structure, ligands, and cofactors
def parse_pdb(pdb_file):
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

    cofactor_names = ["MG", "ZN", "CA", "K", "NA", "FE", "CL", "HOH"]
    ligand_names = ["ACP", "TPS", "TMP", "TPP", "LIG"]

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


#def calculate_distance(atom1, atom2):
#    """Calculate Euclidean distance between two atoms based on their coordinates."""
#    return sum((a - b) ** 2 for a, b in zip(atom1['coord'], atom2['coord'])) ** 0.5


def verify_near_residues(input_pdb, ligand_residue, treshold_distance):
    ligand_atoms = [atom for atom in input_pdb['ligands'] 
                    if (atom['res_name'], atom['res_seq'], atom['chain_id']) == ligand_residue]

    # Combine all atom lists
    all_atoms = input_pdb['chains'] + input_pdb['cofactors'] + input_pdb['ligands']

    # Filter out the ligand atoms from all_atoms
    all_atoms = [atom for atom in all_atoms if (atom['res_name'], atom['res_seq'], atom['chain_id']) != ligand_residue]

    near_residues_dict = []

    for atom in all_atoms:
        # Find the closest ligand atom to the current atom
        min_distance, closest_ligand_atom = min((calculate_distance(atom, ligand_atom), ligand_atom) for ligand_atom in ligand_atoms)

        if min_distance <= treshold_distance:
            info = {
                'molecule_name': atom['res_name'],
                'molecule_number': atom['res_seq'],
                'chain': atom['chain_id'],
                'distance': min_distance,
                'molecule_atom': atom['name'],
                'ligand_atom': closest_ligand_atom['name'],
                'molecule_atom_serial': atom['serial_number'],  # Added the serial number for molecule atom
                'ligand_atom_serial': closest_ligand_atom['serial_number']  # Added the serial number for ligand atom
            }
            near_residues_dict.append(info)

    # Sort the list based on distance
    near_residues_dict.sort(key=lambda x: x['distance'])

    return near_residues_dict


# Define a function to check for potential hydrogen bonds
def is_interaction(atom1_name, atom2_name, residue_name, distance):

    # Check for hydrophobic interactions
    if residue_name in hydrophobic_residues:
        if atom1_name in hydrophobic_atoms or atom2_name in hydrophobic_atoms:
            if distance <= hydrophobic_distance_threshold:
                return "Hydrophobic"
            
        # Specifically check for a carbon-hydrogen interaction
        if (atom1_name.startswith("C") and atom2_name == "H") or (atom2_name.startswith("C") and atom1_name == "H"):
            if distance <= hydrophobic_distance_threshold:
                return "Hydrophobic"
            
    # Check for ionic interaction
    if distance < 4:
        if residue_name in ionic_interactions:
            if atom1_name.startswith(tuple(ionic_interactions.get(residue_name, []))) or \
                atom2_name.startswith(tuple(ionic_interactions.get(residue_name, []))):
                return "Ionic interaction"
            
        elif (atom1_name in ionic_atoms['positive'] and atom2_name in ionic_atoms['negative']) or \
           (atom2_name in ionic_atoms['positive'] and atom1_name in ionic_atoms['negative']):
            return "Ionic interaction"
        
    # Check for hydrogen bond      
    if distance < 3.7:   
        if atom1_name.startswith(tuple(hydrogen_bond_acceptors)) and \
           atom2_name.startswith(tuple(hydrogen_bond_acceptors)) and \
           atom1_name != atom2_name:        
            return "Hydrogen bond"

    # Assuming lennard_jones_potential is defined elsewhere in your code
    v_lj = lennard_jones_potential(atom1_name, atom2_name, residue_name, distance)

    # Check for potential van der Waals interaction based on Lennard-Jones potential
    # Valor Conservador: Se você deseja ser mais conservador e focar apenas nas interações 
    # mais fortes de van der Waals, pode considerar um valor limiar de −0.5 kcal/mol ou mais negativo.
    
    # Valor Moderado: Um valor de −0.2 a −0.3 kcal/mol pode ser uma abordagem intermediária, 
    # onde você identifica interações que têm uma contribuição notável, mas não são extremamente fracas.

    # Análise Detalhada: Se o objetivo é uma análise mais detalhada e abrangente das interações, 
    # incluindo as mais fracas, então −0.1 kcal/mol ou até um pouco mais positivo pode ser aceitável. 
    # No entanto, essas interações devem ser interpretadas com cautela e corroboradas com outras evidências ou análises.

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

    type1 = aminoacid_group_dict[mapped_group1]['PRIMARY MMF TYPE']
    type2 = aminoacid_group_dict[mapped_group2]['PRIMARY MMF TYPE']

    # Get epsilon and sigma values for the atoms
    epsilon1 = float(aminoacid_vdw_dict[type1]['alpha-i'])
    epsilon2 = float(aminoacid_vdw_dict[type2]['alpha-i'])

    sigma1 = float(aminoacid_vdw_dict[type1]['N-i'])
    sigma2 = float(aminoacid_vdw_dict[type2]['N-i'])

    # Combine the epsilon and sigma values
    # Este cálculo refere-se à combinação dos parâmetros de profundidade do poço de energia epsilon para dois átomos 
    # diferentes quando se modela uma interação via potencial de Lennard-Jones. 
    # A combinação geométrica (média geométrica) é comum para este parâmetro.
    epsilon_combined = (epsilon1 * epsilon2) ** 0.5

    # Este cálculo refere-se à combinação dos parâmetros de distância de sigma para os mesmos dois átomos. 
    # Sigma é geralmente interpretado como a distância em que o potencial interatômico entre dois átomos neutros é zero.
    # A combinação aritmética (média aritmética) é típica para este parâmetro.
    sigma_combined = (sigma1 + sigma2) / 2.0

    # Calculate the Lennard-Jones potential
    return 4 * epsilon_combined * ((sigma_combined / r)**12 - (sigma_combined / r)**6)

def map_to_molecular_group(atom_name, residue_name):
    residue_name = residue_name.upper()
    atom_name = atom_name.upper()

    # Verificação de Carbonos
    # Verificação de Carbonos
    if atom_name.startswith("C"):
        carbon_mappings = {
            ("VAL", ("CG1", "CG2")): "CR",
            ("LEU", ("CG", "CD1", "CD2")): "CR",
            ("ILE", ("CG1", "CG2", "CD1")): "CR",
            ("MET", ("CG", "SD", "CE")): "CR",
            ("PHE", ("CZ",)): "CB",
            ("TYR", ("CZ",)): "CB",
            ("TRP", ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")): "CB",
            ("ASP", ("CG", "CD")): "CO2M",
            ("GLU", ("CG", "CD")): "CO2M",
            ("ARG", ("CG",)): "CGD"
        }
        for key, value in carbon_mappings.items():
            if residue_name == key[0] and atom_name in key[1]:
                return value
        if atom_name == "CA" and residue_name != "PRO":
            return "HC"
        
        else: return "CR"
        
    # Verificação de Oxigênios
    if atom_name.startswith("O"):
        if atom_name in ["OD1", "OE1", "OE2", "O"]:
            return "O=C"
        elif atom_name in ["OD2", "OG", "OG1", "OH"]:
            return "OR"
        elif residue_name in ["ASP", "GLU"] and atom_name in ["OD1", "OD2", "OE1", "OE2"]:
            return "O2CM"
        else: 
            return "OR"


    # Verificação de Nitrogênios
    if atom_name.startswith("N"):
        nitrogen_mappings = {
            ("LYS", "NZ"): "NR",
            ("ARG", "NE", "NH1", "NH2"): "NGD+",
            ("HIS", "ND1", "NE2"): "NPYD",
            ("ASN", "ND2"): "NC=O",
            ("GLN", "NE2"): "NC=O",
            ("TRP", "NE1"): "NPYL",
            ("PRO", "N"): "NC=C",
            ("CYS", "NSP"): "NSP"
        }
        for key, value in nitrogen_mappings.items():
            if residue_name == key[0]:
                if atom_name in key[1]:
                    return value
            else:
                return "NR"

    # Verificação de Hidrogênios
    if atom_name.startswith("H"):
        hydrogen_mappings = {
            ("SER", "HO"): "HOR",
            ("THR", "HO"): "HOR",
            ("LYS", "HNZ"): "HNR",
            ("ASP", "HOCO"): "HOCO",
            ("GLU", "HOCO"): "HOCO"
        }
        for key, value in hydrogen_mappings.items():
            if residue_name == key[0]:
                if atom_name in key[1]:
                    return value
        if atom_name == "HN" or atom_name[1:].isdigit():
            return "HNCO"
        else:
            return "H"

    # Verificação de Enxofres
    if atom_name.startswith("S"):
        sulfur_mappings = {
            ("MET", ["SD"]): "S",
            ("CYS", ["SG"]): "HS"
        }
        for key, value in sulfur_mappings.items():
            if residue_name == key[0]:
                if atom_name in key[1]:
                    return value
                
    if atom_name.startswith("MG"):
        return "MG+2"
                
    print(f"O átomo {atom_name} não foi reconhecido.")
    return "UNKNOWN"




def format_line(atom_data, atom_type="ATOM  "):
    return (
        f"{atom_type:6s}{atom_data['serial_number']:5d} {atom_data['name']:<4s} {atom_data['alt_loc']:1s}{atom_data['res_name']:<3s} "
        f"{atom_data['chain_id']:1s}{atom_data['res_seq']:4d}{atom_data['icode']:1s}   "
        f"{atom_data['coord'][0]:8.3f}{atom_data['coord'][1]:8.3f}{atom_data['coord'][2]:8.3f}"
        f"{atom_data['occupancy']:6.2f}{atom_data['temp_factor']:6.2f}          "
        f"{atom_data['element']:^2s}{atom_data['charge']:2s}\n"
    )

def print_pdb_structure(pdb_dict):
    """
    Prints the PDB structure from the dictionary in an organized manner.

    Parameters:
        pdb_dict (dict): Dictionary containing the parsed PDB lists.

    Returns:
        None: Simply prints the PDB structured data.
    """
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

#########################################################################################


# Função para calcular a distância entre dois átomos
def calculate_distance(atom1, atom2):
    coord1 = atom1['coord']
    coord2 = atom2['coord']
    return math.sqrt(sum([(c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)]))

# Função para calcular o produto interno de dois vetores
def dot_product(v1, v2):
    return sum([a*b for a, b in zip(v1, v2)])

# Função para calcular a norma (magnitude) de um vetor
def norm(v):
    return math.sqrt(dot_product(v, v))

# Função para calcular o produto vetorial de dois vetores
def cross_product(v1, v2):
    return [
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    ]

# Função para calcular o ângulo entre três pontos
def calculate_angle(coord_A, coord_B, coord_C):

    # 1 - Obtemos os vetores BA e BC.
    BA = [a-b for a, b in zip(coord_A, coord_B)]
    BC = [c-b for c, b in zip(coord_C, coord_B)]

    # 2 - Calculamos o produto escalar (dot product) desses dois vetores.
    # 3 - Calculamos a magnitude (norma) de cada vetor.
    # 4 - O cosseno do ângulo entre os vetores é dado por: cos(θ) = BA⋅BC / |BA|⋅|BC|
    cosine_angle = dot_product(BA, BC) / (norm(BA) * norm(BC))
    
    # 5 - Usamos a função arco-cosseno (acos) para obter o ângulo θ em radianos.
    angle = math.acos(max(-1.0, min(1.0, cosine_angle)))

    # 6 - Convertemos o ângulo de radianos para graus.
    return math.degrees(angle)

# Função para calcular o diedro entre quatro pontos
def calculate_dihedral(coord_A, coord_B, coord_C, coord_D):
    
    # 1 - Calculamos os vetores BA, CB e DC.
    # 2 - O primeiro plano é definido pelos vetores BA e CB, e o segundo plano é definido pelos vetores CB e DC.

    BA = [a-b for a, b in zip(coord_A, coord_B)]
    CB = [b-c for b, c in zip(coord_B, coord_C)]
    DC = [c-d for c, d in zip(coord_C, coord_D)]

    # 3 - Calculamos os vetores normais a esses planos usando o produto vetorial (cross product). 
    # O vetor normal ao primeiro plano é:  N1 = BA x CB    
    # E o vetor normal ao segundo plano é: N2 = CB x DC.
    normal1 = cross_product(BA, CB)
    normal2 = cross_product(CB, DC)
    n1_norm = norm(normal1)
    n2_norm = norm(normal2)

    # 4 - Calculamos o cosseno do ângulo entre os vetores normais usando o produto escalar e as magnitudes:
    #        cos(θ) = N1⋅N2 / |N1|⋅|N2|
    # Além disso, para garantir que o diedro esteja no intervalo correto de -180° a 180°, 
    # levamos em consideração a direção do vetor formado pelo produto vetorial dos vetores normais e o vetor CB. 
    # Se a direção for negativa, invertemos o sinal do cosseno.
    if n1_norm != 0:
        normal1 = [n/n1_norm for n in normal1]
    if n2_norm != 0:
        normal2 = [n/n2_norm for n in normal2]
    cosine_angle = dot_product(normal1, normal2)
    direction = dot_product(cross_product(normal1, normal2), CB)
    if direction < 0:
        cosine_angle = -cosine_angle

    # 5 - Usamos a função arco-cosseno (acos) para obter o ângulo θ em radianos.
    # 6 - Convertemos o ângulo de radianos para graus.
    angle = math.acos(max(-1.0, min(1.0, cosine_angle)))

    return math.degrees(angle)



# Função para extrair coordenadas de átomos com base em seus números de série
def extract_coordinates(parsed_data, serial_numbers):
    all_atoms = parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands']
    serial_to_coord = {atom['serial_number']: atom['coord'] for atom in all_atoms}
    return [serial_to_coord[serial] for serial in serial_numbers if serial in serial_to_coord]


def extract_atom_info(parsed_data, serial_numbers):
    """
    Retrieve the coordinates and atom names based on their serial numbers.
    """
    all_atoms = parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands']
    serial_to_info = {atom['serial_number']: (atom['coord'], atom['name'] + str(atom['res_seq'])) for atom in all_atoms}
    return [serial_to_info[serial] for serial in serial_numbers if serial in serial_to_info]


def calculate_angles_for_atoms(parsed_data, near_residues_dict):
    
    # Para calcular um ângulo, precisamos de três pontos. O ângulo é formado entre o primeiro e o terceiro ponto, usando o segundo ponto como o vértice.
    # Se tivermos os pontos A, B e C, o ângulo que queremos calcular é o formado por  BA e BC.

    serial_numbers = [entry['molecule_atom_serial'] for entry in near_residues_dict]
    atom_info = extract_atom_info(parsed_data, serial_numbers)
    coordinates = [info[0] for info in atom_info]
    atom_names = [info[1] for info in atom_info]
    
    angles = []
    for i in range(len(coordinates) - 2):
        angle = calculate_angle(coordinates[i], coordinates[i+1], coordinates[i+2])
        angles.append({
            'Atoms': (serial_numbers[i], serial_numbers[i+1], serial_numbers[i+2]),
            'Angle (degrees)': angle,
            'Atom Names': (atom_names[i], atom_names[i+1], atom_names[i+2])
        })
    return angles


def calculate_dihedrals_for_atoms(parsed_data, near_residues_dict):
    
    # Para calcular um diedro, precisamos de quatro pontos. O diedro é o ângulo entre dois planos, e esses planos são definidos pelos pontos.
    # Se tivermos os pontos A, B, C e D, queremos calcular o ângulo entre o plano formado por A, B e C e o plano formado por B, C e D.

    serial_numbers = [entry['molecule_atom_serial'] for entry in near_residues_dict]
    atom_info = extract_atom_info(parsed_data, serial_numbers)
    coordinates = [info[0] for info in atom_info]
    atom_names = [info[1] for info in atom_info]
    
    dihedrals = []
    for i in range(len(coordinates) - 3):
        dihedral = calculate_dihedral(coordinates[i], coordinates[i+1], coordinates[i+2], coordinates[i+3])
        dihedrals.append({
            'Atoms': (serial_numbers[i], serial_numbers[i+1], serial_numbers[i+2], serial_numbers[i+3]),
            'Dihedral (degrees)': dihedral,
            'Atom Names': (atom_names[i], atom_names[i+1], atom_names[i+2], atom_names[i+3])
        })
    return dihedrals



def calculate_angles_for_nearest_atoms(parsed_data, near_residues_dict):
    all_atoms = parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands']

    def find_two_nearest_atoms(target_atom, atoms_list):
        distances = [(atom, calculate_distance(target_atom, atom)) for atom in atoms_list if atom != target_atom]
        distances.sort(key=lambda x: x[1])
        return [distances[0][0], distances[1][0]]

    angles = []
    for entry in near_residues_dict:
        target_atom_serial = entry['molecule_atom_serial']
        target_atom = next(atom for atom in all_atoms if atom['serial_number'] == target_atom_serial)
        
        nearest_atoms = find_two_nearest_atoms(target_atom, all_atoms)
        
        angle = calculate_angle(nearest_atoms[0]['coord'], target_atom['coord'], nearest_atoms[1]['coord'])
        
        angles.append({
            'Target Atom'     : target_atom_serial,
            'Nearest Atoms'   : (nearest_atoms[0]['serial_number'], nearest_atoms[1]['serial_number']),
            'Angle (degrees)' : angle,
            'Atom Names'      : (nearest_atoms[0]['res_name'] + "("+ nearest_atoms[0]['name'] + ")" +  str(nearest_atoms[0]['res_seq']), 
                                target_atom['res_name']       + "("+ target_atom['name']      + ")" +  str(target_atom['res_seq']),
                                nearest_atoms[1]['res_name']  + "("+ nearest_atoms[1]['name'] + ")" +  str(nearest_atoms[1]['res_seq']))
        })


    return angles


def calculate_dihedrals_for_nearest_atoms(parsed_data, near_residues_dict):
    all_atoms = parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands']

    def find_three_nearest_atoms(target_atom, atoms_list):
        distances = [(atom, calculate_distance(target_atom, atom)) for atom in atoms_list if atom != target_atom]
        distances.sort(key=lambda x: x[1])
        return [distances[0][0], distances[1][0], distances[2][0]]

    dihedrals = []
    for entry in near_residues_dict:
        target_atom_serial = entry['molecule_atom_serial']
        target_atom = next(atom for atom in all_atoms if atom['serial_number'] == target_atom_serial)
        
        nearest_atoms = find_three_nearest_atoms(target_atom, all_atoms)
        
        dihedral = calculate_dihedral(nearest_atoms[0]['coord'], nearest_atoms[1]['coord'], nearest_atoms[2]['coord'], target_atom['coord'])
        
        dihedrals.append({
            'Target Atom': target_atom_serial,
            'Nearest Atoms': (nearest_atoms[0]['serial_number'], nearest_atoms[1]['serial_number'], nearest_atoms[2]['serial_number']),
            'Dihedral (degrees)': dihedral,
            'Atom Names': ("("+nearest_atoms[0]['name']+")" + str(nearest_atoms[0]['res_seq']), 
                           "("+nearest_atoms[1]['name']+")" + str(nearest_atoms[1]['res_seq']), 
                           "("+nearest_atoms[2]['name']+")" + str(nearest_atoms[2]['res_seq']),
                           "("+target_atom['name']+")"      + str(target_atom['res_seq']))
        })

    return dihedrals



def get_distance_between_atoms(parsed_data, atom1_idx, atom2_idx):
    atom1 = next(atom for atom in parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands'] if atom['serial_number'] == atom1_idx)
    atom2 = next(atom for atom in parsed_data['chains'] + parsed_data['cofactors'] + parsed_data['ligands'] if atom['serial_number'] == atom2_idx)
    return calculate_distance(atom1, atom2)



def set_output_angle_dihedral(near_residues_dict, ligand_residue_tuple, parsed_data, output_name):
    ligand_name, ligand_num, ligand_chain = ligand_residue_tuple
    
    # Ensure that angles_data and dihedrals_data have at least the same length as near_residues_dict
    angles_data = calculate_angles_for_nearest_atoms(parsed_data, near_residues_dict) + [{} for _ in range(len(near_residues_dict))]
    dihedrals_data = calculate_dihedrals_for_nearest_atoms(parsed_data, near_residues_dict) + [{} for _ in range(len(near_residues_dict))]
    
    interacting_molecules_count = 0

    with open(output_name, 'w', newline='') as file:
        writer = csv.writer(file)
        columns = ["Chain", "Nearby atoms", "Interaction", "Distance (Å)", 
                   "Angle (°)", "Angle distance (Å)", "Angle atoms name", 
                   "Dihedral (°)", "Dihedral distance (Å)", "Dihedral atoms name"]
        writer.writerow(columns)
        print("{:^5} {:^30} {:^20} {:^10} {:^10} {:^15} {:^40} {:^20} {:^5} {:^40}".format(*columns))

        for idx, entry in enumerate(near_residues_dict):
            chain_id = entry['chain']
            distance = entry['distance']
            molecule_atom = entry['molecule_atom']
            ligand_atom = entry['ligand_atom']
            aa_name = entry['molecule_name']
            aa_num = entry['molecule_number']

            atom1_str = f"{molecule_atom}({aa_name}{aa_num})"
            atom2_str = f"{ligand_atom}({ligand_name}{ligand_num})"
            nearby_atoms_str = f"{atom1_str:<12}-{atom2_str:>13}"

            probable_interaction = is_interaction(molecule_atom, ligand_atom, aa_name, distance)
            if probable_interaction not in ["Non-specific"]:
                interacting_molecules_count += 1

            angle_data = angles_data[idx]
            angle, angle_atoms, angle_names = round(angle_data.get('Angle (degrees)', ""), 2), angle_data.get('Nearest Atoms', ""), ", ".join(angle_data.get('Atom Names', ""))
            angle_distance = round(get_distance_between_atoms(parsed_data, *angle_atoms[:2]), 2) if 'Nearest Atoms' in angle_data else ""

            dihedral_data = dihedrals_data[idx]
            dihedral, dihedral_atoms, dihedral_names = round(dihedral_data.get('Dihedral (degrees)', ""), 2), dihedral_data.get('Nearest Atoms', ""), ", ".join(dihedral_data.get('Atom Names', ""))
            dihedral_distance = round(get_distance_between_atoms(parsed_data, *dihedral_atoms[1:3]), 2) if 'Nearest Atoms' in dihedral_data else ""

            writer.writerow([chain_id, nearby_atoms_str, probable_interaction, round(distance, 2), round(angle,2), angle_distance, angle_names, dihedral, dihedral_distance, dihedral_names])
            print("{:^5} {:^30} {:^20} {:^10.2f} {:^10} {:^15} {:^50} {:^15} {:^15} {:^50}".format(chain_id, nearby_atoms_str, probable_interaction, distance, angle, angle_distance, angle_names, dihedral, dihedral_distance, dihedral_names))

    print("\nTotal number of interacting molecules:", interacting_molecules_count)
    print("\n")
    return f"Successfully processed and saved! Total interactions: {interacting_molecules_count}"



if __name__ == "__main__":

    # Arquivos de entrada e saida a serem fornecidos
    input_pdb      = "./3c9t.pdb"  #sys.argv[1]    # EXAMPLE.pdb
    input_molecule = "ACP"         #sys.argv[2]    # ATP
    output_name    = "out.csv"     #sys.argv[3]    # ATP_OUT (csv)
 
    # Distance from the selected molecule
    treshold_distance = 4.0

    # executa e salva o resultados para a classificacao dos contatos
    input_pdb      = parse_pdb(input_pdb)
    ligand_residue = find_molecule(input_pdb, input_molecule)
    near_residues  = verify_near_residues(input_pdb, ligand_residue, treshold_distance)
    set_output_angle_dihedral(near_residues, ligand_residue, input_pdb, output_name)

    # consertar a formatacao da coluna "Dihedral atoms name"