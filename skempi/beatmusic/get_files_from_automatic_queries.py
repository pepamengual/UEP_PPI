from selenium import webdriver
import urllib
import requests
import os.path
import time

def read_file(input_file):
    data = {}
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip().split()
            pdb, chain, link = line[0], line[1], line[2]
            data.setdefault(link, (pdb, chain))
    return data

def runner(link, pdb, chain):
    driver = webdriver.Chrome(executable_path="/Users/pep/beatmusic/chromedriver")
    driver.get(link)
    
    links = driver.find_elements_by_xpath("//a[@href]")
    result = []
    for link in links:
        link = link.get_attribute("href")
        if "tmp/0/" in link:
            result.append(link)
    time.sleep(1)
    driver.stop_client()
    driver.close()
    return result[-1]

def downloader(result):
    r = requests.get(result)
    text = r.text
    text_list = text.split("\n")
    return text_list

def writter(pdb, chain, text_list):
    filename = "output/beatmusic_{}_{}.txt".format(pdb, chain)
    with open(filename, "w") as f:
        for text in text_list:
            f.write(text + "\n")

def main():
    input_file = "chinofarmeo_automatik.txt"
    data = read_file(input_file)
    for link, pdb_chain_tupple in data.items():
        pdb = pdb_chain_tupple[0]
        chain = pdb_chain_tupple[1]
        filename = "output/beatmusic_{}_{}.txt".format(pdb, chain)
        if not os.path.exists(filename):
            print("Working with pdb: {} chain: {}".format(pdb, chain))
            result = runner(link, pdb, chain)
            text_list = downloader(result)
            writter(pdb, chain, text_list)

main()
