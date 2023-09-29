from selenium import webdriver
from selenium.webdriver.common.by import By
import time
import csv

def extract_smiles_from_onclick(onclick_attr):
    # Exemplo de valor: setClipboard('CCC(C)C(NC(=O)C(Cc1ccccc1)NC(=O)C(N)CC(=O)O)C(=O)O')
    prefix = "setClipboard('"
    suffix = "')"
    return onclick_attr[len(prefix):-len(suffix)]

def search_bindingdb(query):
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

    smiles = []

    # Encontra todos os botões com a classe "nonb m1fontSize" e extrai SMILES de seu atributo 'onclick'
    copy_buttons = driver.find_elements(By.XPATH, '//button[contains(@class, "nonb m1fontSize")]')
    for button in copy_buttons:
        onclick_attr = button.get_attribute('onclick')
        if onclick_attr:
            smiles.append(extract_smiles_from_onclick(onclick_attr))
    
    driver.quit()
    
    return smiles

def save_to_csv(data, filename="bindingdb_results.csv"):
    with open(filename, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["SMILES"])  # Write header
        for smile in data:
            writer.writerow([smile])

if __name__ == "__main__":
    query = "N[C@@H](CC(=O)Nc1ccc(Oc2cc(F)c(F)cc2Br)cc1)C(=O)O"
    results = search_bindingdb(query)
    if results:
        save_to_csv(results)
        print(f"Results for {query} have been saved to 'bindingdb_results.csv'.")
    else:
        print(f"No results for {query} found.")
