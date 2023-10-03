from selenium import webdriver
from selenium.webdriver.common.by import By
import time
import csv
import os

def extract_data_from_page(driver):
    """Função para extrair os dados da página e retornar uma lista."""
    data = []
    
    # Encontra todos os botões com a classe "nonb m1fontSize" e extrai SMILES de seu atributo 'onclick'
    copy_buttons = driver.find_elements(By.XPATH, '//button[contains(@class, "nonb m1fontSize")]')
    for button in copy_buttons:
        onclick_attr = button.get_attribute('onclick')
        if onclick_attr and onclick_attr.startswith("setClipboard"):
            # Pega o conteúdo dentro de setClipboard('...') e remove qualquer aspas duplas
            content = onclick_attr.split("'")[1]
            data.append(content)
    
    # Remove entradas que começam com aspas duplas
    data = [item for item in data if not item.startswith('"')]
    return data

import pandas as pd
import os

def save_to_csv(data, filename="bindingdb_results.csv"):
    """Salva os dados no arquivo CSV, ignorando entradas que começam com aspas duplas."""

    # Lê o arquivo existente, se ele existir
    if os.path.exists(filename):
        existing_df = pd.read_csv(filename)
        
        # Remova as linhas que começam com aspas duplas
        existing_df = existing_df[~existing_df.iloc[:, 0].str.startswith('"')]
    else:
        existing_df = pd.DataFrame()

    # Convertendo os dados para DataFrame
    new_df = pd.DataFrame(data)

    # Concatenando os dados existentes com os novos dados
    combined_df = pd.concat([existing_df, new_df], ignore_index=True)

    # Salva tudo em um novo arquivo CSV
    combined_df.to_csv(filename, index=False)





def search_bindingdb(query, similarity_value):
    # Inicializa o WebDriver para Safari
    driver = webdriver.Safari()

    # Vai para a página inicial do BindingDB
    driver.get("https://www.bindingdb.org/rwd/bind/index.jsp")

    # Localiza o campo de pesquisa e insere a query
    input_element = driver.find_element(By.XPATH, '//input[@placeholder="Search by protein (target) name, compound name, author, article title, SMILES, InChi"]')
    input_element.send_keys(query)
    input_element.send_keys(u'\ue007')  # press Enter

    # Aguarda até que os resultados sejam carregados
    time.sleep(10)

    # Encontra o campo de similaridade, limpa, e define o novo valor
    similarity_input = driver.find_element(By.NAME, "Similarity")
    similarity_input.clear()
    similarity_input.send_keys(similarity_value)

    # Clica no botão "GO"
    go_button = driver.find_element(By.XPATH, '//button[text()="GO"]')
    go_button.click()

    time.sleep(10) # Aguarda a página carregar após clicar no botão GO.

    # Extrai os dados
    data = extract_data_from_page(driver)

    # Imprime os dados no prompt
    for item in data:
        print(item)
    
    # Salva os dados em CSV
    save_to_csv(data)
    
    driver.quit()

if __name__ == "__main__":
    query = "N[C@@H](CC(=O)Nc1ccc(Oc2cc(F)c(F)cc2Br)cc1)C(=O)O"
    similarity_value = "0.9"
    search_bindingdb(query, similarity_value)