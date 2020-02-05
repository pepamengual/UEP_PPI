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
    MCC = round(((TP * TN) - (FP * FN)) / ((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))**0.5, 3)

    print("\tP\tN\tPPV/NPV")
    print("P\t{}\t{}\t{}".format(TP, FP, PPV))
    print("N\t{}\t{}\t{}".format(FN, TN, NPV))
    print("\tTPR\tFPR\tMCC")
    print("\t{}\t{}\t{}\n".format(TPR, FPR, MCC))

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
