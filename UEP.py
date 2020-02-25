import argparse
from predictor import read_skempi, scoring_all, compute_statistics, make_models, scoring_single_contact
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
        skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed, skempi_raw_renamed_original, map_beatmusic = read_skempi.process_skempi_data(skempi_path)
        data = make_models.export_mutations(skempi_processed_data_single)
        data_no_renamed = make_models.export_mutations(skempi_processed_data_single_no_renamed)
        
        ### --- MAKE 3D STRUCTURES USING FOLDX --- ###
        #make_models.run_multiprocessing_models(skempi_uep_predictions)
        
        ### ---- UEP ---- ###
        skempi_uep_predictions = scoring_all.run_multiprocessing(skempi_processed_data_single, cpus, training_data)
        uep_results = compute_statistics.mcc(skempi_uep_predictions, skempi_processed_data_single, 1.01, "UEP")
        #compute_statistics.best_mcc(skempi_uep_predictions, skempi_processed_data_single)
        
        training_data_single = load("trained_model/single_contact_matrix", compression="lzma", set_default_extension=False)
        skempi_uep_predictions_single = scoring_single_contact.run_multiprocessing(skempi_processed_data_single, cpus, training_data_single)
        single = {}
        for k, v in skempi_uep_predictions_single.items():
            if k in skempi_uep_predictions:
                single.setdefault(k, v)
        uep_results_2 = compute_statistics.mcc(single, skempi_processed_data_single, 1.01, "UEP_single")

        ### -- PYDOCK --- ###
        #make_models.run_multiprocessing_pydock(data)
        interaction_data_pydock = make_models.get_interaction_data_pydock(data)
        ddG_data_pydock = make_models.get_ddG_pydock(interaction_data_pydock)
        names_pydock = make_models.get_foldx_mutation_names_for_pydock(ddG_data_pydock)
        pydock_results = compute_statistics.mcc(names_pydock, skempi_processed_data_single, 0.0, "pyDock")
        #compute_statistics.best_mcc(names_pydock, skempi_processed_data_single)
        
        ### --- FOLDX --- ###
        #make_models.run_multiprocessing_foldx(data)
        interaction_data_foldx = make_models.get_interaction_data_foldx(data)
        ddG_data_foldx = make_models.get_ddG_foldx(interaction_data_foldx)
        names_foldx = make_models.get_foldx_mutation_names(ddG_data_foldx)
        foldx_results = compute_statistics.mcc(names_foldx, skempi_processed_data_single, 0.0, "FoldX")
        #compute_statistics.best_mcc(names_foldx, skempi_processed_data_single)
        
        ### - CONSENSUS - ###
        consensus = compute_statistics.make_consensus(skempi_uep_predictions, names_pydock, names_foldx, 1.01)
        consensus_results = compute_statistics.mcc(consensus, skempi_processed_data_single, 0.0, "Consensus")

        ### - BEATMUSIC - ###      
        beatmusic_folder = "skempi/beatmusic/output/"
        interaction_data_beatmusic = make_models.run_beatmusic(beatmusic_folder, skempi_uep_predictions, skempi_raw_renamed_original, map_beatmusic)
        beatmusic_results = compute_statistics.mcc(interaction_data_beatmusic, skempi_processed_data_single, 0.0, "BeAtMuSiC")
        #beatmusic_results = compute_statistics.mcc(interaction_data_beatmusic, skempi_processed_data_single_no_renamed, 0.0, "BeAtMuSiC")
        #compute_statistics.best_mcc(interaction_data_beatmusic, skempi_processed_data_single_no_renamed)
        
        ### -- PRODIGY -- ###
        results_prodigy = make_models.run_multiprocessing_prodigy(data)
        prodigy_results = compute_statistics.mcc(results_prodigy, skempi_processed_data_single, 0.0, "PRODIGY")

        ### --- MCSM ---- ###
        mcsm_folder = "skempi/mcsm/output/"
        mcsm_training_path = "skempi/mcsm/dataset/BeAtMuSiC_dataset/BeAtMuSiC.csv"
        interaction_training_data_mcsm, interaction_new_data_mcsm = make_models.run_mcsm_and_split(mcsm_folder, mcsm_training_path, skempi_uep_predictions, skempi_raw_renamed_original)
        mcsm_trained_results = compute_statistics.mcc(interaction_training_data_mcsm, skempi_processed_data_single, 0.0, "mCSM\ntrained")
        mcsm_untrained_results = compute_statistics.mcc(interaction_new_data_mcsm, skempi_processed_data_single, 0.0, "mCSM\nuntrained")
        #compute_statistics.best_mcc(interaction_data_mcsm, skempi_processed_data_single)
        uep_untrained = {}
        for k, v in skempi_uep_predictions.items():
            k = "{}_{}".format(k.split("_")[0], k.split("_")[-1])
            if k in interaction_new_data_mcsm:
                uep_untrained.setdefault(k, v)
        uep_untrained_results = compute_statistics.mcc(uep_untrained, skempi_processed_data_single, 1.01, "UEP\nuntrained")

        all_data = {**uep_results, **pydock_results, **foldx_results, **beatmusic_results, **prodigy_results, **mcsm_trained_results, **mcsm_untrained_results, **consensus_results}
        print(all_data)

        ### -- ALL AGREE --- ###
        compute_statistics.all_agree_matrix(skempi_uep_predictions, names_pydock, names_foldx, results_prodigy, interaction_data_beatmusic, skempi_processed_data_single)
        
        consensus_four = compute_statistics.make_consensus_four(skempi_uep_predictions, names_pydock, names_foldx, results_prodigy, 1.01)
        consensus_results_four = compute_statistics.mcc(consensus_four, skempi_processed_data_single, 0.0, "Consensus")

if __name__ == "__main__":
    cpu, skempi, scan = parse_args()
    main(cpus=cpu, skempi=skempi, scan=scan)
