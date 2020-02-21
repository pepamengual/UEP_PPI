from multiprocessing import Pool
import prody
import itertools
from os import path
import numpy as np

def scoring_skempi(pdb_name, candidate_list, training_data):
    aa_code_singlet_to_triplet = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS',
            'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP',
            'Y':'TYR','V':'VAL'}
    frequency_random_model = get_frequency_random_model(training_data)
    data = {}
    radius = str(7.00)
    path_pdb = "PDBs/{}.pdb".format(pdb_name)
    if not path.exists(path_pdb):
        return None
    pdb = prody.parsePDB(path_pdb)
    for candidates in candidate_list:
        candidate_chunks = candidates.split("_")[-1].split(",")
        counts_original, counts_mutation = 0, 0
        list_to_check = {}
        for candidate in candidate_chunks:
            chain = candidate[1]
            aa_original = candidate[0]
            position = candidate[2:-1]
            aa_mutation = candidate[-1]
            list_to_check.setdefault("{}_{}_{}".format(chain, aa_code_singlet_to_triplet[aa_original], position), aa_code_singlet_to_triplet[aa_mutation])
        for candidate in candidate_chunks:
            chain = candidate[1]
            aa_original = candidate[0]
            position = candidate[2:-1]
            aa_mutation = candidate[-1]
            new_ratio = volume_corr(aa_original, aa_mutation)
            #add changing near residues function for all more than one mutations
            near_residues_original = pdb.select('(ca same residue as within '+ radius +' of (noh resid '+ position +' and chain '+ chain +')) and not chain '+ chain +'')
            near_residues_mutation = pdb.select('(ca same residue as within '+ new_ratio +' of (noh resid '+ position +' and chain '+ chain +')) and not chain '+ chain +'')
            if near_residues_original == None or near_residues_mutation == None:
                 continue

            aa_original_three = (aa_code_singlet_to_triplet[aa_original],)
            aa_mutation_three = (aa_code_singlet_to_triplet[aa_mutation],)
        
            environment_original = tuple(sorted(near_residues_original.getResnames().tolist()))
            if len(candidate_chunks) > 1:
                resnames = near_residues_original.getResnames().tolist()
                chains = near_residues_original.getChids().tolist()
                resnums = near_residues_original.getResnums().tolist()
                for i, r in enumerate(resnames):
                    find = "{}_{}_{}".format(chains[i], resnames[i], resnums[i])
                    if find in list_to_check.keys():
                        resnames[i] = list_to_check[find]

                environment_mutation = tuple(sorted(resnames))
            else:
                environment_mutation = tuple(sorted(near_residues_mutation.getResnames().tolist()))
            if len(environment_original) < 2 or len(environment_mutation) < 2: # or
                continue
            for pair in environment_original:
                if pair in training_data:
                    counts_original += training_data[pair][aa_original_three[0]] * sum(training_data[pair].values())
            for pair in environment_mutation:
                if pair in training_data:
                    counts_mutation += training_data[pair][aa_mutation_three[0]] * sum(training_data[pair].values())
        if counts_mutation != 0 and counts_original != 0:   
            ratio = (counts_mutation/counts_original)# * (frequency_random_model[aa_mutation_three]/frequency_random_model[aa_original_three])
            #if len(candidates.split("_")[1]) == 1 and len(candidates.split("_")[2]) == 1:
            data.setdefault(candidates, round(ratio, 3))
    return data

def volume_corr(aa_original, aa_mutation):
    volume_dict = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140}
    mutation, original = volume_dict[aa_mutation], volume_dict[aa_original]
    min_volume, max_volume = min(volume_dict.values()), max(volume_dict.values())
    diff_volume_min_max = np.abs(max_volume - min_volume)
    diff_volume_mutation_original = np.abs(mutation - original)

    max_adding = 6
    ratio_volume = mutation/original
    volume_modifier = round(diff_volume_mutation_original * max_adding / diff_volume_min_max, 2)
    if ratio_volume > 1:
        new_ratio = str(round(5 + volume_modifier, 2))
    else:
        new_ratio = str(round(5 - volume_modifier, 2))
        if float(new_ratio) < 4:
           new_ratio = str(4)
    return str(7)


def get_frequency_random_model(training_data):
    frequency_random_model = {}
    for environment, target_dict in training_data.items():
        for target, count in target_dict.items():
            frequency_random_model.setdefault(target, 0)
            frequency_random_model[target] += count
    return frequency_random_model

def run_multiprocessing(experimental_skempi_ratios, cpus, training_data):
    skempi_uep_predictions = {}
    pool = Pool(processes=cpus)
    multiple_results = []
    
    experimental_skempi_data = {}
    for candidate, aff_ratio in experimental_skempi_ratios.items():
        pdb = candidate.split("_")[0]
        experimental_skempi_data.setdefault(pdb, {}).setdefault(candidate, aff_ratio)

    for pdb, candidate_dict in experimental_skempi_data.items():
        candidate_list = list(candidate_dict.keys())
        multiple_results.append(pool.apply_async(scoring_skempi, (pdb, candidate_list, training_data)))
    for result in multiple_results:
        data = result.get()
        if data != None:
            for name, prediction in data.items():
                skempi_uep_predictions.setdefault(name, prediction)
    pool.terminate()
    return skempi_uep_predictions

