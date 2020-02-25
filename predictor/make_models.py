import glob
from multiprocessing import Pool
from pathlib import Path
import shutil
import os

def export_mutations(skempi_processed_data_single):
    data = {}
    for key, value in skempi_processed_data_single.items():
        if os.path.exists("PDBs/{}.pdb".format(key.split("_")[0])):
            data.setdefault(key.split("_")[0], []).append(key.split("_")[-1])
    return data

def generate_files(data):
    Path("mutations_foldx").mkdir(parents=True, exist_ok=True)
    for pdb, mutations in data.items():
        if os.path.exists("PDBs/{}.pdb".format(pdb)):
            Path("mutations_foldx/{}".format(pdb)).mkdir(parents=True, exist_ok=True)
            shutil.copy("PDBs/{}.pdb".format(pdb), "mutations_foldx/{}/{}.pdb".format(pdb, pdb))
            shutil.copy("rotabase.txt", "mutations_foldx/{}/rotabase.txt".format(pdb))
            with open("mutations_foldx/{}/individual_list.txt".format(pdb), "w") as f:
                for mutation in mutations:
                    f.write("{};".format(mutation) + "\n")

def make_models_foldx(pdb):
    os.chdir("mutations_foldx/{}".format(pdb))
    os.system("/home/pepamengual/foldx/new_foldx/foldx --command=BuildModel --pdb={}.pdb --mutant-file=individual_list.txt".format(pdb))
    os.chdir("../..")

def interaction_energy(pdb):
    if os.path.exists("mutations_foldx/{}".format(pdb)):
        os.chdir("mutations_foldx/{}".format(pdb))
        all_files = os.listdir()
        for f in all_files:
            if "_" in f and f.endswith(".pdb"):
                os.system("/home/pepamengual/foldx/new_foldx/foldx --command=AnalyseComplex --pdb={}".format(f))
        os.chdir("../..")

def get_interaction_data_foldx(data):
    interaction_data = {}
    pdb_list = list(data.keys())
    for pdb in pdb_list:
        if os.path.exists("mutations_foldx/{}".format(pdb)):
            os.chdir("mutations_foldx/{}".format(pdb))
            all_files = os.listdir()
            for fi in all_files:
                if fi.startswith("Summary"):
                    with open(fi, "r") as f:
                        for line in f:
                            if line.startswith("./"):
                                line = line.rstrip().split()
                                pdb = line[0].split("./")[1]
                                group_1 = line[1]
                                group_2 = line[2]
                                interaction_energy = float(line[5])
                                if pdb.startswith("WT"):
                                    pdb_name = pdb.split("WT_")[1]
                                    interaction_data.setdefault(pdb_name, {}).setdefault("{}_{}".format(group_1, group_2), {}).setdefault("WT", interaction_energy)
                                else:
                                    interaction_data.setdefault(pdb, {}).setdefault("{}_{}".format(group_1, group_2), {}).setdefault("MT", interaction_energy)
            os.chdir("../..")
    return interaction_data

def get_ddG_foldx(interaction_data):
    ddG_data = {}
    for pdb, chain_dict in interaction_data.items():
        ddG_sum = 0
        for chains, wt_mut_dict in chain_dict.items():
            if "WT" in wt_mut_dict and "MT" in wt_mut_dict:
                wt = wt_mut_dict["WT"]
                mt = wt_mut_dict["MT"]
                ddG = mt - wt
                ddG_sum += ddG
        ddG_data.setdefault(pdb, round(-ddG_sum, 5))
    return ddG_data

def get_foldx_mutation_names(ddG_data):
    names = {}
    for name, value in ddG_data.items():
        pdb = name.split("_")[0]
        code = int(name.split("_")[1].split(".pdb")[0])
        with open("mutations_foldx/{}/individual_list.txt".format(pdb), "r") as f:
            for i, line in enumerate(f):
                mutation = line.rstrip().split(";")[0]
                if i + 1 == code:
                    names.setdefault("{}_{}".format(pdb, mutation), value)
    return names

def create_pydock_ini(pdb, mutation_chain, other_chains):
    with open("{}.ini".format(pdb.split(".pdb")[0]), "w") as f:
        f.write("[receptor]" + "\n")
        f.write("pdb     = {}".format(pdb) + "\n")
        f.write("mol     = {}".format(",".join(other_chains)) + "\n")
        f.write("newmol  = {}".format(",".join(other_chains)) + "\n")
        f.write("\n")
        f.write("[ligand]" + "\n")
        f.write("pdb     = {}".format(pdb) + "\n")
        f.write("mol     = {}".format(mutation_chain) + "\n")
        f.write("newmol  = {}".format(mutation_chain) + "\n")

def open_individual_list():
    chain_list = []
    with open("individual_list.txt", "r") as f:
        for line in f:
            mutation = line.rstrip().split(";")[0]
            chain = mutation[1]
            chain_list.append(chain)
    return chain_list

def find_all_chains(pdb):
    all_chains_in_pdb = set()
    with open(pdb, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                all_chains_in_pdb.add(chain)
    return list(all_chains_in_pdb)

def run_pydock(pdb):
    if os.path.exists("mutations_foldx/{}".format(pdb)):
        os.chdir("mutations_foldx/{}".format(pdb))
        all_files = os.listdir()
        chain_list_of_mutations = open_individual_list()
        all_chains_in_pdb = find_all_chains("{}.pdb".format(pdb))
        for f in all_files:
            if "_" in f and f.endswith(".pdb"):
                mutation_chain = chain_list_of_mutations[int(f.split(".pdb")[0].split("_")[-1]) - 1]
                other_chains = list(set(all_chains_in_pdb) - set(mutation_chain))
                create_pydock_ini(f, mutation_chain, other_chains)
                os.system("/home/pepamengual/pydock/pyDock3/pyDock3 {} bindEy".format(f.split(".pdb")[0]))
        os.chdir("../..")

def get_interaction_data_pydock(data):
    interaction_data_pydock = {}
    pdb_list = list(data.keys())
    for pdb in pdb_list:
        if os.path.exists("mutations_foldx/{}".format(pdb)):
            os.chdir("mutations_foldx/{}".format(pdb))
            all_files = os.listdir()
            for fi in all_files:
                if fi.endswith(".ene"):
                    with open(fi, "r") as f:
                        for i, line in enumerate(f):
                            if i == 2:
                                line = line.rstrip().split()
                                pdb = fi
                                interaction_energy = float(line[4])
                                if pdb.startswith("WT"):
                                    pdb_name = pdb.split("WT_")[1]
                                    interaction_data_pydock.setdefault(pdb_name, {}).setdefault("WT", interaction_energy)
                                else:
                                    interaction_data_pydock.setdefault(pdb, {}).setdefault("MT", interaction_energy)
            os.chdir("../..")
    return interaction_data_pydock

def get_ddG_pydock(interaction_data_pydock):
    ddG_data_pydock = {}
    for pdb, wt_mut_dict in interaction_data_pydock.items():
        if "WT" in wt_mut_dict and "MT" in wt_mut_dict:
            wt = wt_mut_dict["WT"]
            mt = wt_mut_dict["MT"]
            ddG = mt - wt
            ddG_data_pydock.setdefault(pdb, round(-ddG, 5))
    return ddG_data_pydock

def get_foldx_mutation_names_for_pydock(ddG_data_pydock):
    names_pydock = {}
    for name, value in ddG_data_pydock.items():
        pdb = name.split("_")[0]
        code = int(name.split("_")[1].split(".ene")[0])
        with open("mutations_foldx/{}/individual_list.txt".format(pdb), "r") as f:
            for i, line in enumerate(f):
                mutation = line.rstrip().split(";")[0]
                if i + 1 == code:
                    names_pydock.setdefault("{}_{}".format(pdb, mutation), value)
    return names_pydock

def run_beatmusic(folder, data_uep, skempi_raw_renamed_original, map_beatmusic): #beatmusic_folder, skempi_uep_predictions, skempi_raw_renamed_original
    #folder = "skempi/beatmusic/output/"
    data_to_test = list(data_uep.keys())
    data_to_test = ["{}_{}".format(i.split("_")[0], i.split("_")[-1]) for i in data_to_test]
    data_to_test = [skempi_raw_renamed_original[i] for i in data_to_test]
    data_to_test = set(data_to_test)
    data = {}
    for path in glob.glob("{}*.txt".format(folder)):
        with open(path, "r") as f:
            for line in f:
                if len(line) > 1:
                    line = line.rstrip().split()
                    pdb, chain, resnum, original, mutation, ddG = line[0].split(".pdb")[0].upper(), path.split("_")[-1].split(".txt")[0], line[3], line[4], line[5], -float(line[7])
                    name = "{}{}{}{}".format(original, chain, resnum, mutation)
                    entry = "{}_{}".format(pdb, name)
                    if len(data_to_test.intersection(set([entry]))) == 1:
                        entry = map_beatmusic[entry] #added this, remove if doesn't work
                        data.setdefault(entry, ddG)
    return data

def run_mcsm(folder, data_uep):
    #folder = "skempi/mcsm/output/"
    data_to_test = list(data_uep.keys())
    data_to_test = ["{}_{}".format(i.split("_")[0], i.split("_")[-1]) for i in data_to_test]
    data_to_test = set(data_to_test)
    data = {}
    for path in glob.glob("{}*.txt".format(folder)):
        with open(path, "r") as f:
            next(f)
            for line in f:
                if len(line) > 1:
                    line = line.rstrip().split()
                    pdb, chain, original, resnum, mutation, ddG = line[0].split(".pdb")[0].upper(), line[1], line[2], line[3], line[4], float(line[6])
                    name = "{}{}{}{}".format(original, chain, resnum, mutation)
                    entry = "{}_{}".format(pdb, name)
                    if len(data_to_test.intersection(set([entry]))) == 1:
                        data.setdefault(entry, ddG)
    return data

def read_mcsm_training_data(training_path):
    #training_path = "skempi/mcsm/dataset/BeAtMuSiC_dataset/BeAtMuSiC.csv"
    mcsm_training_list = []
    with open(training_path, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip().split(";")
            pdb, original, chain, resnum, mutation = line[0].split(".pdb")[0].upper(), line[1], line[2], line[3], line[4] 
            name = "{}_{}{}{}{}".format(pdb, original, chain, resnum, mutation)
            mcsm_training_list.append(name)
    return mcsm_training_list

def run_mcsm_and_split(folder, training_path, data_uep, skempi_raw_renamed_original):
    #data = dict: key = renamed, value = ddG 
    #skempi_raw_original_renamed = dict: key = original, value = renamed
    #data_from_training = list of mutations from training
    data = run_mcsm(folder, data_uep) # folder = "skempi/mcsm/output/"
    mcsm_training_list = read_mcsm_training_data(training_path) # training_path = "skempi/mcsm/dataset/BeAtMuSiC_dataset/BeAtMuSiC.csv"
    training_data, new_data = {}, {}
    for renamed_mutation, ddG in data.items():
        original_mutation = skempi_raw_renamed_original[renamed_mutation]
        if original_mutation in mcsm_training_list:
            training_data.setdefault(renamed_mutation, ddG)
        else:
            new_data.setdefault(renamed_mutation, ddG)
    return training_data, new_data

def run_prodigy(work):
    pdb_mutation_name, file_name_mutation, file_name_wildtype = work[0], work[1], work[2]
    mutation_energy = os.popen("prodigy {} -q".format(file_name_mutation)).read()
    wildtype_energy = os.popen("prodigy {} -q".format(file_name_wildtype)).read()
    if mutation_energy and wildtype_energy:
        wt_energy = float(wildtype_energy.split()[1])
        mt_energy = float(mutation_energy.split()[1])
        if wt_energy and mt_energy:
            ddG = round(wt_energy - mt_energy, 3)
            if ddG is None:
                ddG = "CAN'T"
            return pdb_mutation_name, ddG
    else:
        print("{} can't".format(pdb_mutation_name))
        return pdb_mutation_name, "CAN'T"

def run_multiprocessing_prodigy(data_dict):
    result_file = "skempi/prodigy_results.txt"
    data = {}
    if os.path.isfile(result_file):
         with open(result_file, "r") as f:
            for line in f:
                line = line.rstrip().split()
                mutation, ddG = line[0], float(line[1])
                data.setdefault(mutation, ddG)
    else:
        prodigy_run_list = create_prodigy_jobs(data_dict)
        pool = Pool(processes=25)
        multiple_results = []
        for work in prodigy_run_list:
            multiple_results.append(pool.apply_async(run_prodigy, (work,)))
        for result in multiple_results:
            pdb_mutation_name, ddG = result.get()
            if ddG != "CAN'T":
                data.setdefault(pdb_mutation_name, ddG)

        with open(result_file, "w") as f:
            for pdb_mutation_name, ddG in data.items():
                to_write = "{} {}".format(pdb_mutation_name, ddG)
                f.write(to_write + "\n")
    return data

def run_multiprocessing_models(skempi_processed_data_single):
    data = export_mutations(skempi_processed_data_single)
    generate_files(data)
    pool = Pool(processes=25)
    multiple_results = []
    pdb_list = list(data.keys())
    for pdb in pdb_list:
        multiple_results.append(pool.apply_async(make_models_foldx, (pdb,)))
    for result in multiple_results:
        result.get()
    pool.terminate()

def run_multiprocessing_foldx(data):
    pool = Pool(processes=25)
    multiple_results = []
    pdb_list = list(data.keys())
    for pdb in pdb_list:
        multiple_results.append(pool.apply_async(interaction_energy, (pdb,)))
    for result in multiple_results:
        result.get()
    pool.terminate()

def run_multiprocessing_pydock(data):
    pdb_list = list(data.keys())
    pool = Pool(processes=25)
    multiple_results = []
    for pdb in pdb_list:
        multiple_results.append(pool.apply_async(run_pydock, (pdb,)))
    for result in multiple_results:
        result.get()
    pool.terminate()
