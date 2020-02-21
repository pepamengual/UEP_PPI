from compress_pickle import load

def load_training_data(pickle_path):
    training_data = load(pickle_path, compression="lzma", set_default_extension=False)
    return training_data

def break_down_data(training_data, amino_acid_list, properties):
    data = {"polar": {"polar": 0, "apolar": 0, "+charge": 0, "-charge": 0}, "apolar": {"polar": 0, "apolar": 0, "+charge": 0, "-charge": 0}, "+charge": {"polar": 0, "apolar": 0, "+charge": 0, "-charge": 0}, "-charge": {"polar": 0, "apolar": 0, "+charge": 0, "-charge": 0}}
    for pair, residue_dict in training_data.items():
        res1 = pair[0]
        res2 = pair[1]
        for residue, count in residue_dict.items():
            residue = residue[0]
            if res1 in amino_acid_list and res2 in amino_acid_list and residue in amino_acid_list:
                residue_prop = properties[residue]
                res1_prop = properties[res1]
                res2_prop = properties[res2]
                data[residue_prop][res1_prop] += count
                data[residue_prop][res2_prop] += count
    
    for quality, quality_dict in data.items():
        suma = sum(quality_dict.values())
        for quality_2, value in quality_dict.items():
            fraction = round(value/suma, 2)
            print(quality, quality_2, fraction)

def main():
    amino_acid_list = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", 
                       "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    properties = {"ALA": "apolar", "CYS": "apolar", "ASP": "-charge", "GLU": "-charge", "PHE": "apolar", 
                  "GLY": "apolar", "HIS": "polar", "ILE": "apolar", "LYS": "+charge", "LEU": "apolar", 
                  "MET": "apolar", "ASN": "polar", "PRO": "apolar", "GLN": "polar", "ARG": "+charge", 
                  "SER": "polar", "THR": "polar", "VAL": "apolar", "TRP": "apolar", "TYR": "apolar"}
    pickle_path = "UEP_trained_model_4"
    training_data = load_training_data(pickle_path)
    break_down_data(training_data, amino_acid_list, properties)
main()
