import glob
import os
import prody
from multiprocessing import Pool
from compress_pickle import dump

def training_with_multiprocessing(radius, number_of_processors, path_training_folders):
    training_data = {}
    pool = Pool(processes=number_of_processors)
    multiple_results = []
    for folder in glob.glob(path_training_folders):
        pdb_list = glob.glob(os.path.join(folder, "*.pdb"))
        for pdb_file in pdb_list:
            multiple_results.append(pool.apply_async(environment_creator, (pdb_file, radius)))
    for result in multiple_results:
        observed_contacts_dictionary = result.get()
        for pair_residues, target_residue_dictionary in observed_contacts_dictionary.items():
            if pair_residues in training_data:
                for target_residue, observed_contacts in target_residue_dictionary.items():
                    if target_residue in training_data[pair_residues]:
                        training_data[pair_residues][target_residue] += observed_contacts
                    else:
                        training_data[pair_residues][target_residue] = observed_contacts
            else:
                 training_data[pair_residues] = target_residue_dictionary
    pool.terminate()

    return training_data

def environment_creator(pdb_file, radius):
    observed_contacts_dictionary = {}
    radius = str(radius)
    pdb = prody.parsePDB(pdb_file)
    pdb_chains_list = pdb.getChids()

    for chain in set(pdb_chains_list):
        chain_interface = pdb.select('(ca same residue as within '+ radius +'.00 of (noh chain '+ chain +')) and not chain '+ chain +'')
        if chain_interface is None:
            continue

        for residue_id, number_id, chain_id, icode_id, atom_id in zip(chain_interface.getResnames(), chain_interface.getResnums(), chain_interface.getChids(), chain_interface.getIcodes(), chain_interface.getIndices().tolist()):
            if not icode_id:
                icode_id = "_"
            number_id = '`{}{}`'.format(str(number_id), icode_id)
            first_residue_picker = pdb.select('ca chain '+ chain_id +' and resid '+ number_id +'')
            first_residue_residue_id = list(first_residue_picker.getResnames())[0]

            residue_interface = pdb.select('(ca same residue as within '+ radius +'.00 of (noh resid '+ number_id +' and chain '+ chain_id +')) and not chain '+ chain_id +'')

            if residue_interface is None:
                continue

            residue_id_contact_list = list(residue_interface.getResnames().tolist())
            for res in residue_id_contact_list:
                observed_contacts_dictionary.setdefault(res, {}).setdefault(first_residue_residue_id, 0)
                observed_contacts_dictionary[res][first_residue_residue_id] += 1
    #print(observed_contacts_dictionary)
    return observed_contacts_dictionary

def main():
    path_training_folders="/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*"
    radius = 4
    number_of_processors = 27
    training_data = training_with_multiprocessing(radius, number_of_processors, path_training_folders)
    dump(training_data, "single_contact_matrix", compression="lzma", set_default_extension=False)
main()
