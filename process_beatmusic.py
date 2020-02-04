import glob

def read_file(folder):
    data = {}
    for path in glob.glob("{}*.txt".format(folder)):
        with open(path, "r") as f:
            for line in f:
                if len(line) > 1:
                    line = line.rstrip().split()
                    pdb, chain, resnum, original, mutation, ddG = line[0].split(".pdb")[0].upper(), path.split("_")[-1].split(".txt")[0], line[3], line[4], line[5], float(line[7])
                    name = "{}{}{}{}".format(original, chain, resnum, mutation)
                    data.setdefault(pdb, {}).setdefault(name, ddG)
    return data

def main():
    folder = "skempi/beatmusic/output/"
    data = read_file(folder)
    print(data)
main()
