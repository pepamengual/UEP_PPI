def compute_volume_difference(dictionary):
    volume = {"G": 0.00279, "A": 0.05702, "S": 0.09204, "C": 0.14907,
              "T": 0.19341, "D": 0.21051, "P": 0.22790, "N": 0.22972,
              "V": 0.25674, "E": 0.32837, "Q": 0.34861, "I": 0.37671,
              "H": 0.37694, "L": 0.37876, "M": 0.38872, "K": 0.45363,
              "F": 0.55298, "R": 0.58946, "Y": 0.61150, "W": 0.79351}

    volume_data = {}
    for mutation, prediction in dictionary.items():
        wild_type = mutation.split("_")[1][0]
        mutant = mutation.split("_")[1][-1]
        volume_diff = round(volume[mutant] - volume[wild_type], 3)
        
        if volume_diff > -0.10 and volume_diff < 0.10:
            volume_data.setdefault("Neutral", {}).setdefault(mutation, prediction)
        if volume_diff <= -0.10:
            volume_data.setdefault("Loss", {}).setdefault(mutation, prediction)
        if volume_diff >= 0.10:
            volume_data.setdefault("Gain", {}).setdefault(mutation, prediction)
    return volume_data


def compute_hydrophobicity_difference(dictionary):
    hydrophobicity = {"G": 0.48, "A": 0.62, "S": -0.18, "C": 0.29,
                      "T": -0.05, "D": -0.90, "P": 0.12, "N": -0.78,
                      "V": 1.08, "E": -0.74, "Q": -0.85, "I": 1.38,
                      "H": -0.40, "L": 1.06, "M": 0.64, "K": -1.50,
                      "F": 1.19, "R": -2.53, "Y": 0.26, "W": 0.81}

    hydrophobicity_data = {}
    for mutation, prediction in dictionary.items():
        wild_type = mutation.split("_")[1][0]
        mutant = mutation.split("_")[1][-1]
        hydrophobicity_diff = round(hydrophobicity[mutant] - hydrophobicity[wild_type], 3)

        if hydrophobicity_diff > -0.30 and hydrophobicity_diff < 0.30:
            hydrophobicity_data.setdefault("Neutral", {}).setdefault(mutation, prediction)
        if hydrophobicity_diff <= -0.30:
            hydrophobicity_data.setdefault("Loss", {}).setdefault(mutation, prediction)
        if hydrophobicity_diff >= 0.30:
            hydrophobicity_data.setdefault("Gain", {}).setdefault(mutation, prediction)

    return hydrophobicity_data
