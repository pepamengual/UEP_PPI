import glob
from predictor import read_skempi, scoring_all, compute_statistics, make_models
from compress_pickle import load

def run_beatmusic(folder="skempi/beatmusic/output/"):
    data = {}
    for path in glob.glob("{}*.txt".format(folder)):
        with open(path, "r") as f:
            for line in f:
                if len(line) > 1:
                    line = line.rstrip().split()
                    pdb, chain, resnum, original, mutation, ddG = line[0].split(".pdb")[0].upper(), path.split("_")[-1].split(".txt")[0], line[3], line[4], line[5], -float(line[7])
                    name = "{}{}{}{}".format(original, chain, resnum, mutation)
                    entry = "{}_{}".format(pdb, name)
                    data.setdefault(entry, ddG)
    return beatmusic_data

def main():
    skempi_path = "skempi/skempi_v2.csv"
    folder = "skempi/beatmusic/output/"
    interaction_data_beatmusic = run_beatmusic(folder)


    model_trained = "trained_model/UEP_trained_model_4"
    training_data = load(model_trained, compression="lzma", set_default_extension=False)
    skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed = read_skempi.process_skempi_data(skempi_path)
    
    compute_statistics.mcc(interaction_data_beatmusic, skempi_processed_data_single_no_renamed, 0.0, "BEATMUSIC")
    compute_statistics.best_mcc(interaction_data_beatmusic, skempi_processed_data_single_no_renamed)



main()
