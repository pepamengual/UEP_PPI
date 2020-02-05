import glob
from predictor import read_skempi, scoring_all, compute_statistics, make_models
from compress_pickle import load

def main():
    skempi_path = "skempi/skempi_v2.csv"

    model_trained = "trained_model/UEP_trained_model_4"
    training_data = load(model_trained, compression="lzma", set_default_extension=False)
    skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed = read_skempi.process_skempi_data(skempi_path)
    
    data = {}
    for entry, value in skempi_processed_data_single.items():
        pdb = entry.split("_")[0]
        mutation_info = entry.split("_")[-1]
        chain = mutation_info[1]
        mutation_cleaned = "{}{}".format(mutation_info[0], mutation_info[2:])
        data.setdefault(pdb, []).append([chain, mutation_cleaned])
        
    for pdb, info_list in data.items():
        with open("skempi/mcsm/mutation_lists/{}.txt".format(pdb), "w") as f:
            for mutation in info_list:
                chain = mutation[0]
                mutation_cleaned = mutation[1]
                to_write = "{} {}".format(chain, mutation_cleaned)
                f.write(to_write + "\n")
main()
