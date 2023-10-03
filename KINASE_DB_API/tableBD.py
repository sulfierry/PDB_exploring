from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
import os
import time

def initiate_driver():
    driver = webdriver.Safari()
    return driver

def search_bindingdb(query, driver):
    driver.get("https://www.bindingdb.org/rwd/bind/index.jsp")
    search_box = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, 'input[placeholder="Search by protein (target) name, compound name, author, article title, SMILES, InChi"]'))
    )
    search_box.send_keys(query)
    search_box.send_keys(Keys.ENTER)

def extract_table_with_selenium(driver):
    time.sleep(15)  # Wait for the results to load
    
    try:
        # Find the div containing the table
        div_container = WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.XPATH, "//div[.//table[@class='index_table']]"))
        )
        
        # Now, find the table inside this div
        table = div_container.find_element(By.CLASS_NAME, "index_table")
        
        rows = table.find_elements(By.TAG_NAME, "tr")
        
        header = [th.text for th in rows[0].find_elements(By.TAG_NAME, "th")]
        data = []
        
        for row in rows[1:]:
            data.append([td.text for td in row.find_elements(By.TAG_NAME, "td")])
        
        return pd.DataFrame(data, columns=header)

    except Exception as e:
        print(e)
        return None


def main(query="N[C@@H](CC(=O)Nc1ccc(Oc2cc(F)c(F)cc2Br)cc1)C(=O)O"):
    driver = initiate_driver()
    search_bindingdb(query, driver)
    
    # Extract table using Selenium
    results = extract_table_with_selenium(driver)

    # Save the table as CSV
    if results is not None and not results.empty:
        results.to_csv("bindingdb_table_results.csv", index=False)
        print(f"Table results for {query} have been saved to 'bindingdb_table_results.csv'.")
    else:
        print(f"No table results for {query} found.")

    driver.close()

if __name__ == "__main__":
    main()
