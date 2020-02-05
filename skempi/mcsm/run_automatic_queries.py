import os
from selenium import webdriver
import time
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options

def runner(browser_path, website, pdb_path, mutations_path, submit_path, download_path, pdb, txt):
    ### OPENING THE BROWSER ###
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    driver = webdriver.Chrome(executable_path=browser_path, options=chrome_options)
    driver.get(website)

    ### WAIT, UPLOAD FILES AND CLICK SUBMIT ###
    wait = WebDriverWait(driver, 40)
    wait.until(EC.element_to_be_clickable((By.XPATH, pdb_path)))
    wait.until(EC.element_to_be_clickable((By.XPATH, mutations_path)))
    driver.find_element_by_xpath(pdb_path).send_keys(pdb)
    driver.find_element_by_xpath(mutations_path).send_keys(txt)
    wait.until(EC.element_to_be_clickable((By.XPATH, submit_path)))
    submit_button = driver.find_element_by_xpath(submit_path)
    submit_button.click()

    ### WAIT, CLICK DOWNLOAD AND GET RESULT TEXT ###
    wait.until(EC.element_to_be_clickable((By.XPATH, download_path)))
    download_button = driver.find_element_by_xpath(download_path)
    download_button.click()
    text_simulations = driver.find_element_by_xpath(".//pre")
    results = text_simulations.text
    driver.stop_client()
    driver.close()
    return results

def show_and_save_results(pdb_name, export_path, results):
    ### PRINT THE RESULTS ON THE TERMINAL ###
    print("--> Mutations on {}".format(pdb_name))
    print(results)
    print("--> Exporting results to {}".format(export_path))
    ### EXPORT RESULTS INTO A FILE ###
    with open(export_path, "w") as f:
        f.write(results)
    
def main():
    browser_path = "/Users/pep/UEP_PPI/skempi/mcsm/chromedriver"
    website = "http://biosig.unimelb.edu.au/mcsm/protein_protein"

    pdb_path = "/html/body/div[2]/div[2]/div[3]/div/form/input[1]"
    mutations_path = "/html/body/div[2]/div[2]/div[3]/div/form/input[2]"
    submit_path = "/html/body/div[2]/div[2]/div[3]/div/form/button"
    download_path = "/html/body/div[2]/div[2]/div/div[2]/div[2]/a"
    
    for run_file in os.listdir("mutation_lists/"):
        pdb_name = run_file.split(".txt")[0]
        pdb = "/Users/pep/UEP_PPI/skempi/mcsm/PDBs/{}.pdb".format(pdb_name)
        txt = "/Users/pep/UEP_PPI/skempi/mcsm/mutation_lists/{}.txt".format(pdb_name)
        export_path = "output/mcsm_results_{}.txt".format(pdb_name)
        if not os.path.exists(export_path):
            results = runner(browser_path, website, pdb_path, mutations_path, submit_path, download_path, pdb, txt)
            show_and_save_results(pdb_name, export_path, results)
main()
