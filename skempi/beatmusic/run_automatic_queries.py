from selenium import webdriver
import time
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select

def process_data(input_data):
    data = {}
    with open(input_data, "r") as f:
        for line in f:
            pdb = line.rstrip().split(",")[0]
            chain = line.rstrip().split(",")[1]
            data.setdefault(pdb, []).append(chain)
    return data

def runner(pdb, target_chain):
    driver = webdriver.Chrome(executable_path="/Users/pep/beatmusic/chromedriver")
    website = "http://babylone.ulb.ac.be/beatmusic/query.php"
    driver.get(website)
    
    ### FIRST STEP ###
    pdb_box = driver.find_element_by_name("pdb")
    pdb_box.send_keys(pdb)
    submit_button_1 = driver.find_element_by_css_selector('[value="Submit"]')
    submit_button_1.click()
    
    ### SECOND STEP ###
    wait = WebDriverWait(driver, 100)
    wait.until(EC.presence_of_element_located((By.CLASS_NAME, 'pdbinfo')))
    pdb_info = driver.find_elements_by_css_selector(".pdbinfo")
    pdb_info = pdb_info[0].text
    chain_order = []
    for chain_info in pdb_info.split("\n"):
        if chain_info.startswith("Chain"):
            chain_id = chain_info.split(":")[0].split("Chain ")[1]
            chain_order.append(chain_id)
    
    for i, chain in enumerate(chain_order):
        name = "ch{}".format(i)
        if chain == target_chain:
            value = 1
        else:
            value = 2
        action = "input[type='radio'][name='{}'][value='{}']".format(name, value)
        chain_button = driver.find_element_by_css_selector(action)
        chain_button.click()
    
    submit_button_2 = driver.find_element_by_css_selector('[value="Submit"]')
    submit_button_2.click()
    
    ### THIRD STEP ####
    wait.until(EC.presence_of_element_located((By.CLASS_NAME, 'pdbinfo')))
    pdb_info = driver.find_elements_by_css_selector(".pdbinfo")
    pdb_info = pdb_info[0].text
    select_chain = ""
    for chain_info in pdb_info.split("\n"):
        if chain_info.startswith("Sequence-unique entity"):
            chain_id = chain_info.split()[-1]
            if target_chain in chain_id:
                select_chain += chain_id
    
    select = Select(driver.find_element_by_name('sysmode'))
    select.select_by_value(select_chain)
        
    submit_button_3 = driver.find_element_by_css_selector('[value="Submit"]')
    submit_button_3.click()

    ### FOURTH STEP ###
    result = []
    links = driver.find_elements_by_xpath("//a[@href]")
    for link in links:
        link = link.get_attribute("href")
        if "results.php" in link:
            result.append(link)
    driver.stop_client()
    driver.close()
    return result[-1]

def main():
    input_data = "data.txt"
    data = process_data(input_data)
    
    results = {}
    for pdb, list_of_chains in data.items():
        for target_chain in list_of_chains:
            result = runner(pdb, target_chain)
            print(pdb, target_chain, result)
            results.setdefault(pdb, {}).setdefault(target_chain, result)

    with open("chinofarmeo_automatik.txt", "w") as f:
        for pdb, chains in results.items():
            for chain, result in chains.items():
                to_write = "{} {} {}".format(pdb, chain, result)
                f.write(to_write + "\n")
main()
