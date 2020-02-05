import argparse
from predictor import read_skempi, scoring_all, compute_statistics, make_models
from compress_pickle import load

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

def main(cpus=27, skempi=False, scan=""):
    skempi_path = "skempi/skempi_v2.csv"
    model_trained = "trained_model/UEP_trained_model_4"
    if skempi and scan == "":
        ### --- FILTERING SKEMPI --- ###
        training_data = load(model_trained, compression="lzma", set_default_extension=False)
        skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed, skempi_raw_renamed_original = read_skempi.process_skempi_data(skempi_path)
        data = make_models.export_mutations(skempi_processed_data_single)
        data_no_renamed = make_models.export_mutations(skempi_processed_data_single_no_renamed)
            
        ### --- MAKE 3D STRUCTURES USING FOLDX --- ###
        #make_models.run_multiprocessing_models(skempi_uep_predictions)
        
        ### ---- UEP ---- ###
        #skempi_uep_predictions = scoring_all.run_multiprocessing(skempi_processed_data_single, cpus, training_data)
        #compute_statistics.mcc(skempi_uep_predictions, skempi_processed_data_single, 1.01, "UEP")
        #compute_statistics.best_mcc(skempi_uep_predictions, skempi_processed_data_single)

        ### -- PYDOCK --- ###
        #make_models.run_multiprocessing_pydock(data)
        #interaction_data_pydock = make_models.get_interaction_data_pydock(data)
        #ddG_data_pydock = make_models.get_ddG_pydock(interaction_data_pydock)
        #names_pydock = make_models.get_foldx_mutation_names_for_pydock(ddG_data_pydock)
        #compute_statistics.mcc(names_pydock, skempi_processed_data_single, 0.0, "PYDOCK")
        #compute_statistics.best_mcc(names_pydock, skempi_processed_data_single)
        
        ### --- FOLDX --- ###
        #make_models.run_multiprocessing_foldx(data)
        #interaction_data_foldx = make_models.get_interaction_data_foldx(data)
        #ddG_data_foldx = make_models.get_ddG_foldx(interaction_data_foldx)
        #names_foldx = make_models.get_foldx_mutation_names(ddG_data_foldx)
        #compute_statistics.mcc(names_foldx, skempi_processed_data_single, 0.0, "FOLDX")
        #compute_statistics.best_mcc(names_foldx, skempi_processed_data_single)
        
        ### - BEATMUSIC - ###      
        beatmusic_folder = "skempi/beatmusic/output/"
        interaction_data_beatmusic = make_models.run_beatmusic(beatmusic_folder)
        compute_statistics.mcc(interaction_data_beatmusic, skempi_processed_data_single_no_renamed, 0.0, "BEATMUSIC")
        #compute_statistics.best_mcc(interaction_data_beatmusic, skempi_processed_data_single_no_renamed)
        
        ### --- MCSM ---- ###
        mcsm_folder = "skempi/mcsm/output/"
        mcsm_training_path = "skempi/mcsm/dataset/BeAtMuSiC_dataset/BeAtMuSiC.csv"
        interaction_training_data_mcsm, interaction_new_data_mcsm = make_models.run_mcsm_and_split(mcsm_folder, mcsm_training_path, skempi_raw_renamed_original)
        compute_statistics.mcc(interaction_training_data_mcsm, skempi_processed_data_single, 0.0, "TRANING MCSM")
        compute_statistics.mcc(interaction_new_data_mcsm, skempi_processed_data_single, 0.0, "NEW MCSM")
        #compute_statistics.best_mcc(interaction_data_mcsm, skempi_processed_data_single)
if __name__ == "__main__":
    cpu, skempi, scan = parse_args()
    main(cpus=cpu, skempi=skempi, scan=scan)
