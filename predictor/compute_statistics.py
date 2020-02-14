import numpy as np

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
    
    PPV = round(TP/(TP + FP), 3)
    NPV = round(TN/(FN + TN), 3)
    TPR = round(TP/(TP + FN), 3)
    FPR = round(TN/(FP + TN), 3)
    MCC = round(((TP * TN) - (FP * FN)) / ((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))**0.5, 2)

    print("\tP\tN\tPPV/NPV")
    print("P\t{}\t{}\t{}".format(TP, FP, PPV))
    print("N\t{}\t{}\t{}".format(FN, TN, NPV))
    print("\tTPR\tFPR\tMCC")
    print("\t{}\t{}\t{}\n".format(TPR, FPR, MCC))
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

def make_consensus(uep_results, pydock_results, foldx_results):
    consensus_results = {}
    for mutation, uep_value in uep_results.items():
        mutation = "{}_{}".format(mutation.split("_")[0], mutation.split("_")[-1])
        if mutation in pydock_results and mutation in foldx_results:
            pydock_value = pydock_results[mutation]
            foldx_value = foldx_results[mutation]
            positive_count = len([i for i in [uep_value-1.01, pydock_value, foldx_value] if i > 0])
            if positive_count >= 2:
                consensus_results.setdefault(mutation, 2)
            else:
                consensus_results.setdefault(mutation, -0.5)
    return consensus_results

def correlations(uep_results, pydock_results, foldx_results, beatmusic_results, prodigy_results):
    correlation_data = {"UEP-pyDock": 0, "UEP-FoldX": 0, "pyDock-FoldX": 0, "UEP-BeAtMuSiC": 0, "UEP-PRODIGY": 0, "UEP-pyDock-FoldX": 0}
    for mutation, uep_value in uep_results.items():
        mutation = "{}_{}".format(mutation.split("_")[0], mutation.split("_")[-1])
        if mutation in pydock_results and mutation in foldx_results and mutation in prodigy_results:
            uep_value -= 1.01 # to compare it with the others, UEP ratio is at 1
            pydock_value = pydock_results[mutation]
            foldx_value = foldx_results[mutation]
            #beatmusic_value = beatmusic_results[mutation]
            prodigy_value = prodigy_results[mutation]
            if (uep_value > 0 and pydock_value > 0) or (uep_value <= 0 and pydock_value <= 0):
                correlation_data["UEP-pyDock"] += 1
            if (uep_value > 0 and foldx_value > 0) or (uep_value <= 0 and foldx_value <= 0):
                correlation_data["UEP-FoldX"] += 1
            if (pydock_value > 0 and foldx_value > 0) or (pydock_value <= 0 and foldx_value <= 0):
                correlation_data["pyDock-FoldX"] += 1
            #if (uep_value > 0 and beatmusic_value > 0) or (uep_value <= 0 and beatmusic_value <= 0):
            #    correlation_data["UEP-BeAtMuSiC"] += 1
            if (uep_value > 0 and prodigy_value > 0) or (uep_value <= 0 and prodigy_value <= 0):
                correlation_data["UEP-PRODIGY"] += 1
            if (uep_value > 0 and pydock_value > 0 and foldx_value > 0) or (uep_value <= 0 and pydock_value <= 0 and foldx_value <= 0):
                correlation_data["UEP-pyDock-FoldX"] += 1
    print(correlation_data)

