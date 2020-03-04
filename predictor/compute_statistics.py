import numpy as np
from functools import reduce

def mcc(skempi_uep_predictions, experimental_skempi_ratios, threshold, name):
    print("------> {}".format(name))
    experimental_threshold = 0
    confusion_matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    
    experimental_skempi_ratios_shorter = {}
    for candidate, score in experimental_skempi_ratios.items():
        #if len(candidate.split("_")[1]) > 1 and len(candidate.split("_")[2]) > 1:
        #    continue
        candidate = "{}_{}".format(candidate.split("_")[0], candidate.split("_")[-1])
        experimental_skempi_ratios_shorter.setdefault(candidate, score)
    
    for candidate, uep_score in skempi_uep_predictions.items():
        candidate = "{}_{}".format(candidate.split("_")[0], candidate.split("_")[-1])
        if candidate in experimental_skempi_ratios_shorter:
            experimental_ratio = experimental_skempi_ratios_shorter[candidate]
        
            if uep_score > threshold and experimental_ratio > experimental_threshold:
                confusion_matrix["TP"] += 1
            if uep_score > threshold and experimental_ratio <= - experimental_threshold:
                confusion_matrix["FP"] += 1
            if uep_score <= threshold and experimental_ratio > experimental_threshold:
                confusion_matrix["FN"] += 1
            if uep_score <= threshold and experimental_ratio <= - experimental_threshold:
                confusion_matrix["TN"] += 1

    TP = confusion_matrix["TP"]
    FP = confusion_matrix["FP"]
    FN = confusion_matrix["FN"]
    TN = confusion_matrix["TN"]
    
    PPV = round(TP/(TP + FP), 2)
    NPV = round(TN/(FN + TN), 2)
    TPR = round(TP/(TP + FN), 2)
    FPR = round(TN/(FP + TN), 2)
    MCC = round(((TP * TN) - (FP * FN)) / ((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))**0.5, 2)

    print("\tP\tN\tPPV/NPV")
    print("P\t{}\t{}\t{}".format(TP, FP, PPV))
    print("N\t{}\t{}\t{}".format(FN, TN, NPV))
    print("\tTPR\tFPR\tMCC\tcounts")
    print("\t{}\t{}\t{}\t{}\n".format(TPR, FPR, MCC, TP+TN+FP+FN))
    return {name: [TP, FP, MCC, FN, TN]}

def best_mcc(skempi_uep_predictions, experimental_skempi_ratios):
    uep_thresholds = np.arange(-2, 2, 0.05)
    mcc = 0
    best_threshold = 0
    experimental_threshold = 0
    
    experimental_skempi_ratios_shorter = {}
    for candidate, score in experimental_skempi_ratios.items():
        #if len(candidate.split("_")[1]) > 1 and len(candidate.split("_")[2]) > 1:
        #    continue
        candidate = "{}_{}".format(candidate.split("_")[0], candidate.split("_")[-1])
        experimental_skempi_ratios_shorter.setdefault(candidate, score)

    for uep_threshold in uep_thresholds:
        confusion_matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
        for candidate, uep_score in skempi_uep_predictions.items():
            candidate = "{}_{}".format(candidate.split("_")[0], candidate.split("_")[-1])
            if candidate in experimental_skempi_ratios_shorter:
                experimental_ratio = experimental_skempi_ratios_shorter[candidate]
                if uep_score > uep_threshold and experimental_ratio > experimental_threshold:
                    confusion_matrix["TP"] += 1
                if uep_score > uep_threshold and experimental_ratio <= - experimental_threshold:
                    confusion_matrix["FP"] += 1
                if uep_score <= uep_threshold and experimental_ratio > experimental_threshold:
                    confusion_matrix["FN"] += 1
                if uep_score <= uep_threshold and experimental_ratio <= - experimental_threshold:
                    confusion_matrix["TN"] += 1
        TP = confusion_matrix["TP"]
        FP = confusion_matrix["FP"]
        FN = confusion_matrix["FN"]
        TN = confusion_matrix["TN"]
        if TP == 0 or FP == 0 or FN == 0 or TN == 0:
            MCC = 0
        else:
            MCC = round(((TP * TN) - (FP * FN)) / ((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))**0.5, 4)
        #print(MCC, uep_threshold)
        if MCC > mcc:
            mcc = MCC
            best_threshold = uep_threshold
    print("Best MCC: {}, Threshold: {}".format(mcc, round(best_threshold, 3)))

def make_consensus(uep_results, pydock_results, foldx_results, value):
    consensus_results = {}
    for mutation, uep_value in uep_results.items():
        mutation = "{}_{}".format(mutation.split("_")[0], mutation.split("_")[-1])
        if mutation in pydock_results and mutation in foldx_results:
            pydock_value = pydock_results[mutation]
            foldx_value = foldx_results[mutation]
            positive_count = len([i for i in [uep_value-value, pydock_value, foldx_value] if i > 0])
            if positive_count >= 2:
                consensus_results.setdefault(mutation, 2)
            else:
                consensus_results.setdefault(mutation, -0.5)
    return consensus_results

def make_consensus_four(uep_results, pydock_results, foldx_results, prodigy_results, value):
    consensus_results = {}
    for mutation, uep_value in uep_results.items():
        mutation = "{}_{}".format(mutation.split("_")[0], mutation.split("_")[-1])
        if mutation in pydock_results and mutation in foldx_results and mutation in prodigy_results:
            pydock_value = pydock_results[mutation]
            foldx_value = foldx_results[mutation]
            prodigy_value = prodigy_results[mutation]
            positive_count = len([i for i in [uep_value-value, pydock_value, foldx_value, prodigy_value] if i > 0])
            if positive_count >= 3:
                consensus_results.setdefault(mutation, 2)
            else:
                consensus_results.setdefault(mutation, -0.5)
    return consensus_results

def sub_all_agree(obligatory_predictor, other_predictors, experimental, predictions_classified):
    experimental_improving = list(experimental["improving"])
    experimental_decreasing = list(experimental["decreasing"])
    obligatory_predicted_improving = list(predictions_classified[obligatory_predictor]["improving"])
    obligatory_predicted_decreasing = list(predictions_classified[obligatory_predictor]["decreasing"])


    tp_lists = [experimental_improving, obligatory_predicted_improving]
    tn_lists = [experimental_decreasing, obligatory_predicted_decreasing]
    fp_lists = [experimental_decreasing, obligatory_predicted_improving]
    fn_lists = [experimental_improving, obligatory_predicted_decreasing]

    for other_predictor in other_predictors:
        predicted_improving = list(predictions_classified[other_predictor]["improving"])
        predicted_decreasing = list(predictions_classified[other_predictor]["decreasing"])
        tp_lists.append(predicted_improving)
        tn_lists.append(predicted_decreasing)
        fp_lists.append(predicted_improving)
        fn_lists.append(predicted_decreasing)

    tp_all = reduce(set.intersection, [set(l_) for l_ in tp_lists])
    tn_all = reduce(set.intersection, [set(l_) for l_ in tn_lists])
    fp_all = reduce(set.intersection, [set(l_) for l_ in fp_lists])
    fn_all = reduce(set.intersection, [set(l_) for l_ in fn_lists])
    
    TP, TN, FP, FN = len(tp_all), len(tn_all), len(fp_all), len(fn_all)
    MCC = round(((TP * TN) - (FP * FN)) / ((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))**0.5, 3)

    PPV, NPV, TPR, FPR = round(TP/(TP + FP), 3), round(TN/(FN + TN), 3), round(TP/(TP + FN), 3), round(TN/(FP + TN), 3)
    
    print("{},".format(obligatory_predictor), ", ".join(other_predictors))
    print("\tP\tN\tPPV/NPV")
    print("P\t{}\t{}\t{}".format(TP, FP, PPV))
    print("N\t{}\t{}\t{}".format(FN, TN, NPV))
    print("\tTPR\tFPR\tMCC\tcounts")
    print("\t{}\t{}\t{}\t{}\n".format(TPR, FPR, MCC, TP+TN+FP+FN))

def all_agree_matrix(uep_results, pydock_results, foldx_results, prodigy_results, beatmusic_results, experimental_skempi_ratios):
    predictions = {"UEP": uep_results, "pyDock": pydock_results, "FoldX": foldx_results, "PRODIGY": prodigy_results, "BeAtMuSiC": beatmusic_results}
    experimental = {"improving": set(), "decreasing": set()}
    for candidate, score in experimental_skempi_ratios.items():
        candidate = "{}_{}".format(candidate.split("_")[0], candidate.split("_")[-1])
        if score > 0:
            experimental["improving"].add(candidate)
        else:
            experimental["decreasing"].add(candidate)

    predictions_classified = {}
    for predictor, prediction_dict in predictions.items():
        for mutation, score in prediction_dict.items():
            mutation = "{}_{}".format(mutation.split("_")[0], mutation.split("_")[-1])
            if predictor == "UEP":
                score -= 1.01
            if score > 0:
                predictions_classified.setdefault(predictor, {}).setdefault("improving", set()).add(mutation)
            else:
                predictions_classified.setdefault(predictor, {}).setdefault("decreasing", set()).add(mutation)

    
    sub_all_agree("PRODIGY", ["FoldX"], experimental, predictions_classified)
    sub_all_agree("pyDock", ["PRODIGY"], experimental, predictions_classified)
    sub_all_agree("pyDock", ["FoldX"], experimental, predictions_classified)
    sub_all_agree("UEP", ["PRODIGY"], experimental, predictions_classified)
    sub_all_agree("UEP", ["pyDock"], experimental, predictions_classified)
    sub_all_agree("UEP", ["FoldX"], experimental, predictions_classified)
    sub_all_agree("pyDock", ["FoldX", "PRODIGY"], experimental, predictions_classified)
    sub_all_agree("UEP", ["FoldX", "PRODIGY"], experimental, predictions_classified)
    sub_all_agree("UEP", ["pyDock", "PRODIGY"], experimental, predictions_classified)
    sub_all_agree("UEP", ["pyDock", "FoldX"], experimental, predictions_classified)
    sub_all_agree("UEP", ["pyDock", "FoldX", "PRODIGY"], experimental, predictions_classified)
    sub_all_agree("UEP", ["BeAtMuSiC"], experimental, predictions_classified)
    sub_all_agree("PRODIGY", ["BeAtMuSiC"], experimental, predictions_classified)
    sub_all_agree("FoldX", ["BeAtMuSiC"], experimental, predictions_classified)
    sub_all_agree("pyDock", ["BeAtMuSiC"], experimental, predictions_classified)
