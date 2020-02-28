from predictor import read_skempi, scoring_without_normalization, compute_statistics, make_models, scoring_single_contact, volume_classifier
from compress_pickle import load

def classify_mutations(data_dictionary):
    data = volume_classifier.compute_volume_difference(data_dictionary)
    gain = {**data["Large_gain"], **data["Medium_gain"]}
    loss = {**data["Large_loss"], **data["Medium_loss"]}
    neutral = {**data["Neutral"], **data["Small_gain"], **data["Small_loss"]}
    return gain, loss, neutral

def score_mutations(gain, neutral, loss, threshold, name, skempi_processed_data_single):
    compute_statistics.mcc(gain, skempi_processed_data_single, threshold, "Gain_{}".format(name))
    compute_statistics.mcc(neutral, skempi_processed_data_single, threshold, "Neutral_{}".format(name))
    compute_statistics.mcc(loss, skempi_processed_data_single, threshold, "Loss_{}".format(name))

def main():
    skempi_path = "skempi/skempi_v2.csv"
    model_trained = "trained_model/UEP_trained_model_4"
    training_data = load(model_trained, compression="lzma", set_default_extension=False)
    skempi_processed_data_single, skempi_processed_data_multiple, skempi_processed_data_single_no_renamed, skempi_raw_renamed_original, map_beatmusic = read_skempi.process_skempi_data(skempi_path)
    ### ---- UEP ---- ###
    skempi_uep_predictions = scoring_without_normalization.run_multiprocessing(skempi_processed_data_single, 27, training_data)
    uep_results = compute_statistics.mcc(skempi_uep_predictions, skempi_processed_data_single, 1.01, "UEP")
    
    mutation_list = {}
    for mutation, predicted_value in skempi_uep_predictions.items():
        if mutation in skempi_processed_data_single:
            mutation = "{}_{}".format(mutation.split("_")[0], mutation.split("_")[-1])
            mutation_list.setdefault(mutation, predicted_value)
    gain, loss, neutral = classify_mutations(mutation_list)
    score_mutations(gain, neutral, loss, 1.01, "UEP", skempi_processed_data_single)

    data = make_models.export_mutations(skempi_processed_data_single)
    interaction_data_pydock = make_models.get_interaction_data_pydock(data)
    ddG_data_pydock = make_models.get_ddG_pydock(interaction_data_pydock)
    names_pydock = make_models.get_foldx_mutation_names_for_pydock(ddG_data_pydock)
    pydock_results = compute_statistics.mcc(names_pydock, skempi_processed_data_single, 0.0, "pyDock")
    gain, loss, neutral = classify_mutations(names_pydock)
    score_mutations(gain, neutral, loss, 0, "pyDock", skempi_processed_data_single)

    data = make_models.export_mutations(skempi_processed_data_single)
    interaction_data_foldx = make_models.get_interaction_data_foldx(data)
    ddG_data_foldx = make_models.get_ddG_foldx(interaction_data_foldx)
    names_foldx = make_models.get_foldx_mutation_names(ddG_data_foldx)
    foldx_results = compute_statistics.mcc(names_foldx, skempi_processed_data_single, 0.0, "FoldX")
    gain, loss, neutral = classify_mutations(names_foldx)
    score_mutations(gain, neutral, loss, 0, "FoldX", skempi_processed_data_single)

    data = make_models.export_mutations(skempi_processed_data_single)
    results_prodigy = make_models.run_multiprocessing_prodigy(data)
    prodigy_results = compute_statistics.mcc(results_prodigy, skempi_processed_data_single, 0.0, "PRODIGY")
    gain, loss, neutral = classify_mutations(results_prodigy)
    score_mutations(gain, neutral, loss, 0, "PRODIGY", skempi_processed_data_single)

    beatmusic_folder = "skempi/beatmusic/output/"
    interaction_data_beatmusic = make_models.run_beatmusic(beatmusic_folder, skempi_uep_predictions, skempi_raw_renamed_original, map_beatmusic)
    beatmusic_results = compute_statistics.mcc(interaction_data_beatmusic, skempi_processed_data_single, 0.0, "BeAtMuSiC")
    gain, loss, neutral = classify_mutations(interaction_data_beatmusic)
    score_mutations(gain, neutral, loss, 0, "BeAtMuSiC", skempi_processed_data_single)

main()   
