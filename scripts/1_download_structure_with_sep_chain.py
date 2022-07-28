""" Import libraries """
from Bio.PDB import *

import pandas as pd

import requests
import os
import time
import shutil
import argparse

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")
log_path = os.path.join(current_path, "logs")

""" Set auxiliary functions """
def download_PDB(structure_id):
	response = requests.get(f"https://files.rcsb.org/download/{structure_id}.pdb")
	return response.text

def process_pdb_file(path_to_read, chain, path_to_save):
	class ChainSelect(Select):
		def __init__(self, needed_chain):
			self.needed_chain = needed_chain

		def accept_chain(self, chain):
			if chain.id == self.needed_chain:
				return 1
			else:
				return 0

	structure_id = path_to_read.split('.')[0]

	parser = PDBParser()
	io = PDBIO()
	structure = parser.get_structure(structure_id, path_to_read)
	if len(structure) > 1:
		input("CHECK")
        # Raise Error #
	else:
		io.set_structure(structure)
		io.save(path_to_save, ChainSelect(chain))

""" Set main function """
def main():
	parser = argparse.ArgumentParser(description='### Here should be the description! ###')
	parser.add_argument("-input", help="The name of your PDB ID with separate chain, for example '4GAW_A'", type=str)
	args = parser.parse_args()
	input_name = args.input

	uncorrected_structure_files = []
	#input_name = input("Enter the name of {PDB ID}_{chain}:\n")
	structure, chain = input_name.split('_')[0], input_name.split('_')[1]

	""" Download structure file """
	structure_data = download_PDB(structure)
	if "ATOM" not in structure_data:
		uncorrected_structure_files.append(f"{structure}.pdb is uncorrected!")
	else:
		list_to_write = []
		for stroka in structure_data.strip('\n').split('\n'):
			list_to_write.append(stroka)
			if stroka.split(' ')[0] == "ENDMDL":
				break
		with open(os.path.join(script_path, f"{structure}.pdb"), "w") as file:
			file.write('\n'.join(list_to_write))

	""" Processing structure file if it is from PDB database """
	process_pdb_file(os.path.join(script_path, f"{structure}.pdb"), chain, os.path.join(script_path, f"{input_name}.pdb"))

	""" Write structures uncorrected """
	if len(uncorrected_structure_files) > 0:
		if not os.path.exists(log_path):
			os.mkdir(log_path)
		with open(os.path.join(log_path, "uncorrected_structure.txt"), "w") as file:
			file.write('\n'.join(uncorrected_structure_files))

""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
