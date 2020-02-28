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

        if abs(volume_diff) < 0.05:
            volume_data.setdefault("Neutral", {}).setdefault(mutation, prediction)
        
        if volume_diff <= -0.05 and volume_diff > -0.10:
            volume_data.setdefault("Small_loss", {}).setdefault(mutation, prediction)
        if volume_diff >= 0.05 and volume_diff < 0.10:
            volume_data.setdefault("Small_gain", {}).setdefault(mutation, prediction)
        
        if volume_diff <= -0.10 and volume_diff > -0.20:
            volume_data.setdefault("Medium_loss", {}).setdefault(mutation, prediction)
        if volume_diff >= 0.10 and volume_diff < 0.20:
            volume_data.setdefault("Medium_gain", {}).setdefault(mutation, prediction)
        
        if volume_diff <= -0.20:
            volume_data.setdefault("Large_loss", {}).setdefault(mutation, prediction)
        if volume_diff >= 0.20:
            volume_data.setdefault("Large_gain", {}).setdefault(mutation, prediction)
    
    neutral = len(volume_data["Neutral"])
    large_loss, medium_loss, small_loss = len(volume_data["Large_loss"]), len(volume_data["Medium_loss"]), len(volume_data["Small_loss"])
    large_gain, medium_gain, small_gain = len(volume_data["Large_gain"]), len(volume_data["Medium_gain"]), len(volume_data["Small_gain"])
    total = neutral + large_loss + medium_loss + small_loss + neutral + small_gain + medium_gain + large_gain

    print(large_loss, medium_loss, small_loss, neutral, small_gain, medium_gain, large_gain)
    print(int(100*large_loss/total), int(100*medium_loss/total), int(100*small_loss/total), int(100*neutral/total), int(100*small_gain/total), int(100*medium_gain/total), int(100*large_gain/total))
    return volume_data
