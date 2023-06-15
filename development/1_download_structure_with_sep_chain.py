""" Import libraries """
from Bio.PDB import *

import pandas as pd

import requests
import os
import sys
import glob
import time
import shutil

''' Set paths '''
main_path = "/data2/Protease/TEST"
input_path = os.path.join(main_path, "input")
blastp_path = "/data2/Protease/TEST/input/blastp"

dataset_path = os.path.join(main_path, "datasets")
if not os.path.exists(dataset_path):
    os.mkdir(dataset_path)

name_of_set = input(f"Enter the name of structure set: {', '.join([i for i in os.listdir(blastp_path) if '.' not in i])}\n")
StructureSet_path = os.path.join(dataset_path, name_of_set)
if not os.path.exists(StructureSet_path):
    os.mkdir(StructureSet_path)

''' Initialise input files '''
input_file = name_of_set.split('_')[0] + '_structures_' + name_of_set.split('_')[1] + '.csv'
input_data = pd.read_csv(os.path.join(input_path, input_file))
structure_list = input_data["structure"].tolist()

def download_PDB(structure_id):
	response = requests.get(f"https://files.rcsb.org/download/{structure_id}.pdb")
	return response.text

def download_AF(structure_id):
	response = requests.get(f"https://alphafold.ebi.ac.uk/files/{structure_id}-model_v1.pdb")
	return response.text

def process_pdb_file(path_to_read, structure, path_to_save):
    class ChainSelect(Select):

        def __init__(self, needed_chain):
            self.needed_chain = needed_chain

        def accept_chain(self, chain):
            if chain.id == self.needed_chain:
                return 1
            else:
                return 0

    structure_id = structure.split('_')[0]
    structure_chain = structure.split('_')[1]

    parser = PDBParser()
    io = PDBIO()
    structure_io = parser.get_structure(structure_id, f"{path_to_read}/{structure_id}.pdb")
    if len(structure_io) > 1:
        input("CHECK!")
        sys.exit()
    else:
        io.set_structure(structure_io)
        io.save(f"{path_to_save}/{structure}.pdb", ChainSelect(structure_chain))

def main():
    uncorrected_structure_files = []
    for num, structure in enumerate(structure_list):
        print(f"{num + 1} --- {structure}")
        structure_id = structure.split('_')[0] if '_' in structure else structure
        
        structure_path = os.path.join(StructureSet_path, structure_id)
        if not os.path.exists(structure_path):
            os.mkdir(structure_path)

        ''' Download structure file '''
        structure_file = os.path.join(structure_path, f"{structure_id}.pdb")
        if not os.path.exists(structure_file):
            structure_data = download_PDB(structure_id) if "AF-" not in structure_id else download_AF(structure_id)
            if "ATOM" not in structure_data:
                uncorrected_structure_files.append(f"{structure}")
                shutil.rmtree(structure_path)
                continue
            else:
                list_to_write = []
                for stroka in structure_data.strip('\n').split('\n'):
                    list_to_write.append(stroka)
                    if stroka.split(' ')[0] == "ENDMDL":
                        break
                with open(structure_file, "w") as file:
                    file.write('\n'.join(list_to_write))

        ''' Processing structure file if it is from PDB database '''
        if "AF-" not in structure_id:
            process_pdb_file(structure_path, structure, structure_path)

    ''' Write structures uncorrected '''
    if len(uncorrected_structure_files) > 0:
        with open(f"{StructureSet_path}/{name_of_set}_UncorrectedStructureFiles.txt", "w") as file:
            file.write('\n'.join(uncorrected_structure_files))

''' Launch script '''
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
