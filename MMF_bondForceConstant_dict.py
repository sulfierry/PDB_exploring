#                Copyright (c) Merck and Co., Inc., 1994
#                         All Rights Reserved

# Transcribed by: Leon Sulfierry GitHub:https://github.com/sulfierry

"""
Neste novo dicionário, que se refere às constantes de força de ligação no MMFF94, as chaves e os valores associados têm os seguintes significados:

**Chave (tupla de dois valores)**:
1. **Primeiro valor**: Representa o tipo de átomo em uma extremidade da ligação.
2. **Segundo valor**: Representa o tipo de átomo na outra extremidade da ligação.

Então, a chave identifica uma ligação específica entre dois tipos de átomos no modelo de força de campo. Por exemplo, a chave ('1', '6') representa a ligação entre átomos do tipo 1 e tipo 6.

**Valores associados à chave**:
1. **'r0-ref'**: Este é o comprimento de ligação de referência (ou equilíbrio) para os dois tipos de átomos. É a distância ideal entre os átomos quando a energia potencial da ligação é mínima. As unidades são tipicamente em Ångströms (Å).
2. **'kb-ref'**: É a constante de força para a ligação, que determina a rigidez da ligação. Esta constante é usada para calcular a energia potencial associada a desvios do comprimento de ligação de equilíbrio. As unidades são geralmente em kcal/mol·Å².
3. **'Source'**: Indica a fonte ou referência de onde esses valores foram obtidos ou ajustados. 'C94' e 'E94' são provavelmente referências ou códigos para estudos ou publicações específicas que forneceram ou validaram esses parâmetros.
"""

bndk_dict = {
    ('1', '6'): {'r0-ref': '1.084', 'kb-ref': '5.15', 'Source': 'C94'},
    ('1', '7'): {'r0-ref': '1.001', 'kb-ref': '7.35', 'Source': 'C94'},
    ('1', '8'): {'r0-ref': '0.947', 'kb-ref': '9.10', 'Source': 'C94'},
    ('1', '9'): {'r0-ref': '0.92', 'kb-ref': '10.6', 'Source': 'E94'},
    ('1', '14'): {'r0-ref': '1.48', 'kb-ref': '2.3', 'Source': 'E94'},
    ('1', '15'): {'r0-ref': '1.415', 'kb-ref': '2.95', 'Source': 'C94'},
    ('1', '16'): {'r0-ref': '1.326', 'kb-ref': '4.30', 'Source': 'C94'},
    ('1', '17'): {'r0-ref': '1.28', 'kb-ref': '4.3', 'Source': 'E94'},
    ('1', '35'): {'r0-ref': '1.41', 'kb-ref': '4.2', 'Source': 'E94'},
    ('1', '53'): {'r0-ref': '1.60', 'kb-ref': '2.7', 'Source': 'E94'},
    ('6', '6'): {'r0-ref': '1.512', 'kb-ref': '3.80', 'Source': 'C94'},
    ('6', '7'): {'r0-ref': '1.439', 'kb-ref': '4.55', 'Source': 'C94'},
    ('6', '8'): {'r0-ref': '1.393', 'kb-ref': '5.40', 'Source': 'C94'},
    ('6', '9'): {'r0-ref': '1.353', 'kb-ref': '6.20', 'Source': 'C94'},
    ('6', '14'): {'r0-ref': '1.86', 'kb-ref': '2.6', 'Source': 'E94'},
    ('6', '15'): {'r0-ref': '1.84', 'kb-ref': '2.7', 'Source': 'E94'},
    ('6', '16'): {'r0-ref': '1.812', 'kb-ref': '2.85', 'Source': 'C94'},
    ('6', '17'): {'r0-ref': '1.781', 'kb-ref': '2.75', 'Source': 'C94'},
    ('6', '35'): {'r0-ref': '1.94', 'kb-ref': '2.6', 'Source': 'E94'},
    ('6', '53'): {'r0-ref': '2.16', 'kb-ref': '1.4', 'Source': 'E94'},
    ('7', '7'): {'r0-ref': '1.283', 'kb-ref': '6.00', 'Source': 'C94'},
    ('7', '8'): {'r0-ref': '1.333', 'kb-ref': '5.90', 'Source': 'C94'},
    ('7', '9'): {'r0-ref': '1.36', 'kb-ref': '5.9', 'Source': 'E94'},
    ('7', '14'): {'r0-ref': '1.74', 'kb-ref': '3.7', 'Source': 'E94'},
    ('7', '15'): {'r0-ref': '1.65', 'kb-ref': '4.8', 'Source': 'E94'},
    ('7', '16'): {'r0-ref': '1.674', 'kb-ref': '3.75', 'Source': 'C94'},
    ('7', '17'): {'r0-ref': '1.75', 'kb-ref': '3.5', 'Source': 'E94'},
    ('7', '35'): {'r0-ref': '1.90', 'kb-ref': '2.9', 'Source': 'E94'},
    ('7', '53'): {'r0-ref': '2.10', 'kb-ref': '1.6', 'Source': 'E94'},
    ('8', '8'): {'r0-ref': '1.48', 'kb-ref': '3.6', 'Source': 'E94'},
    ('8', '9'): {'r0-ref': '1.42', 'kb-ref': '4.6', 'Source': 'E94'},
    ('8', '14'): {'r0-ref': '1.63', 'kb-ref': '5.2', 'Source': 'E94'},
    ('8', '15'): {'r0-ref': '1.66', 'kb-ref': '4.7', 'Source': 'E94'},
    ('8', '16'): {'r0-ref': '1.470', 'kb-ref': '9.90', 'Source': 'C94'},
    ('8', '17'): {'r0-ref': '1.70', 'kb-ref': '4.1', 'Source': 'E94'},
    ('8', '35'): {'r0-ref': '1.85', 'kb-ref': '3.4', 'Source': 'E94'},
    ('8', '53'): {'r0-ref': '2.05', 'kb-ref': '1.6', 'Source': 'E94'},
    ('9', '14'): {'r0-ref': '1.57', 'kb-ref': '6.4', 'Source': 'E94'},
    ('9', '15'): {'r0-ref': '1.54', 'kb-ref': '7.1', 'Source': 'E94'},
    ('9', '16'): {'r0-ref': '1.55', 'kb-ref': '6.9', 'Source': 'E94'},
    ('14', '14'): {'r0-ref': '2.32', 'kb-ref': '1.3', 'Source': 'E94'},
    ('14', '15'): {'r0-ref': '2.25', 'kb-ref': '1.5', 'Source': 'E94'},
    ('14', '16'): {'r0-ref': '2.15', 'kb-ref': '2.0', 'Source': 'E94'},
    ('14', '17'): {'r0-ref': '2.02', 'kb-ref': '3.1', 'Source': 'E94'},
    ('14', '35'): {'r0-ref': '2.19', 'kb-ref': '2.1', 'Source': 'E94'},
    ('14', '53'): {'r0-ref': '2.44', 'kb-ref': '1.5', 'Source': 'E94'},
    ('15', '15'): {'r0-ref': '2.21', 'kb-ref': '1.7', 'Source': 'E94'},
    ('15', '16'): {'r0-ref': '2.10', 'kb-ref': '2.4', 'Source': 'E94'},
    ('15', '17'): {'r0-ref': '2.03', 'kb-ref': '3.0', 'Source': 'E94'},
    ('15', '35'): {'r0-ref': '2.21', 'kb-ref': '2.0', 'Source': 'E94'},
    ('15', '53'): {'r0-ref': '2.47', 'kb-ref': '1.4', 'Source': 'E94'},
    ('16', '16'): {'r0-ref': '2.052', 'kb-ref': '2.50', 'Source': 'C94'},
    ('16', '17'): {'r0-ref': '2.04', 'kb-ref': '2.9', 'Source': 'E94'},
    ('16', '35'): {'r0-ref': '2.24', 'kb-ref': '1.9', 'Source': 'E94'},
    ('16', '53'): {'r0-ref': '2.40', 'kb-ref': '1.7', 'Source': 'E94'},
    ('17', '17'): {'r0-ref': '1.99', 'kb-ref': '3.5', 'Source': 'E94'},
    ('35', '35'): {'r0-ref': '2.28', 'kb-ref': '2.4', 'Source': 'E94'},
    ('53', '53'): {'r0-ref': '2.67', 'kb-ref': '1.6', 'Source': 'E94'},
}
