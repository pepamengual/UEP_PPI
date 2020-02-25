import numpy as np

def read_file(input_file):
    skempi_raw_data_single, skempi_raw_data_multiple, skempi_raw_data_single_no_renamed, skempi_raw_renamed_original, map_beatmusic = {}, {}, {}, {}, {}
    with open(input_file, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip().split(";")
            pdb, mutation_not_cleaned, mutation_cleaned = line[0], line[1], line[2]
            id_mutation, id_mutation_not_cleaned = "{}_{}".format(pdb, mutation_cleaned), "{}_{}".format(pdb, mutation_not_cleaned)
            aff_mut_parsed, aff_wt_parsed = line[7], line[9]
            if (aff_mut_parsed and aff_wt_parsed):
                aff_ratio = np.log(float(aff_wt_parsed)/float(aff_mut_parsed))
                aff_ratio = round(aff_ratio, 2)
                if ("," not in id_mutation and id_mutation[-1] != "A"):
                    skempi_raw_data_single.setdefault(id_mutation, []).append(aff_ratio)
                    skempi_raw_data_single_no_renamed.setdefault(id_mutation_not_cleaned, []).append(aff_ratio)
                    map_beatmusic.setdefault("{}_{}".format(id_mutation_not_cleaned.split("_")[0], id_mutation_not_cleaned.split("_")[-1]), "{}_{}".format(id_mutation.split("_")[0], id_mutation.split("_")[-1]))
                    skempi_raw_renamed_original.setdefault("{}_{}".format(id_mutation.split("_")[0], id_mutation.split("_")[-1]), "{}_{}".format(id_mutation_not_cleaned.split("_")[0], id_mutation_not_cleaned.split("_")[-1]))
                if ("," in id_mutation):
                    new_mutations = [i[-1] for i in id_mutation.split(",")]
                    if not "A" in new_mutations:
                        skempi_raw_data_multiple.setdefault(id_mutation, []).append(aff_ratio)
    for k, v in sorted(map_beatmusic.items()):
        print(k, v)

    return skempi_raw_data_single, skempi_raw_data_multiple, skempi_raw_data_single_no_renamed, skempi_raw_renamed_original, map_beatmusic

def filter_redundancy(dictionary):
    dictionary_filtered = {}
    for id_mutation, aff_ratio_list in dictionary.items():
        if all(i > 0 for i in aff_ratio_list) or all(i < 0 for i in aff_ratio_list):
            aff_ratio_mean = round(np.mean(aff_ratio_list), 2)
            dictionary_filtered.setdefault(id_mutation, aff_ratio_mean)
    return dictionary_filtered

def process_skempi_data(input_file):
    skempi_raw_data_single, skempi_raw_data_multiple, skempi_raw_data_single_no_renamed, skempi_raw_renamed_original, map_beatmusic = read_file(input_file)
    skempi_processed_data_single = filter_redundancy(skempi_raw_data_single)
    skempi_processed_data_multiple = filter_redundancy(skempi_raw_data_multiple)
    skempi_processed_data_single_no_renamed = filter_redundancy(skempi_raw_data_single_no_renamed)

    return skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed, skempi_raw_renamed_original, map_beatmusic
