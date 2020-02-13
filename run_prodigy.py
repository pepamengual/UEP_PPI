from multiprocessing import Pool
import argparse
from predictor import read_skempi, scoring_all, compute_statistics, make_models
from compress_pickle import load
import os

HELP = " \
Command:\n \
----------\n \
run skempi benchmark: python3 UEP.py --cpu 27 --skempi True \n \
scan interface of PDB file: python3 UEP.py --cpu 27 --scan PDB.pdb \
"

def parse_args():
    parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu", type=int, help="Amount of CPUs to use", default=27)
    parser.add_argument("--skempi", help="Re-run skempi results", action="store_true")
    parser.add_argument("--scan", type=str, help="Scan a PDB file and generate .csv", default="")
    args = parser.parse_args()
    return args.cpu, args.skempi, args.scan

def create_prodigy_jobs(data):
    work_list = []
    pdb_list = list(set([name.split("_")[0] for name in data.keys()]))
    for pdb in pdb_list:
        mutation_list_path = "mutations_foldx/{}/individual_list.txt".format(pdb)
        if os.path.exists(mutation_list_path):
            with open("mutations_foldx/{}/individual_list.txt".format(pdb)) as f:
                for i, mutation in enumerate(f, start=1):
                    mutation = mutation.split(";")[0]
                    file_name_mutation = "mutations_foldx/{}/{}_{}.pdb".format(pdb, pdb, i)
                    file_name_wildtype = "mutations_foldx/{}/WT_{}_{}.pdb".format(pdb, pdb, i)
                    if os.path.exists(file_name_mutation) and os.path.exists(file_name_wildtype):
                        pdb_mutation_name = "{}_{}".format(pdb, mutation)
                        work_list.append((pdb_mutation_name, file_name_mutation, file_name_wildtype))
    return work_list

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
    prodigy_run_list = create_prodigy_jobs(data_dict)
    data = {}
    pool = Pool(processes=25)
    multiple_results = []
    for work in prodigy_run_list:
        multiple_results.append(pool.apply_async(run_prodigy, (work,)))
    for result in multiple_results:
        pdb_mutation_name, ddG = result.get()
        if ddG != "CAN'T":
            print(pdb_mutation_name, ddG)
            data.setdefault(pdb_mutation_name, ddG)
    
    with open("skempi/prodigy_results.txt", "w") as f:
        for pdb_mutation_name, ddG in data.items():
            to_write = "{} {}".format(pdb_mutation_name, ddG)
            f.write(to_write + "\n")
    return data

def main(cpus=27, skempi=False, scan=""):
    skempi_path = "skempi/skempi_v2.csv"
    model_trained = "trained_model/UEP_trained_model_4"
    if skempi and scan == "":
        ### --- FILTERING SKEMPI --- ###
        training_data = load(model_trained, compression="lzma", set_default_extension=False)
        skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed, skempi_raw_renamed_original = read_skempi.process_skempi_data(skempi_path)
        data = make_models.export_mutations(skempi_processed_data_single)
        results_prodigy = run_multiprocessing_prodigy(data)
        compute_statistics.mcc(results_prodigy, skempi_processed_data_single, 0.0, "PRODIGY")

if __name__ == "__main__":
    cpu, skempi, scan = parse_args()
    main(cpus=cpu, skempi=skempi, scan=scan)
