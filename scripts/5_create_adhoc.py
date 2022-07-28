""" Import libraries """
import os
import time
from itertools import chain as ch

import pandas as pd
import numpy as np

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")
log_path = os.path.join(current_path, "logs")

""" Set important variables """
SS = ['a','b','c','d','e','f','g','h',
	'i','j','k','l','m','n','o','p','q',
	'r','s','t','u','v','w','z','x','y']
HTR = 'X'
b_sheet = ['E']
a_helix = ['H', 'G', 'I']
loop = ['T', 'S', 'B']

""" Set auxiliary functions """
def find_ACC_bp(stroka):
    for index,element in enumerate(stroka):
        if ',' in element:
            break
    return stroka[index-3],stroka[index-2],stroka[index-1]

""" Set auxiliary function """

def form_dataset(path_to_read, name_file_dssp, path_to_save):
	structure_chain = name_file_dssp.split('.')[0]
	structure = structure_chain.split('_')[0] if '_' in structure_chain else structure_chain

	with open(f"{path_to_read}/{name_file_dssp}", 'r') as file:
		data = file.read()

	start_index = [i for i,j in enumerate(data.strip('\n').split('\n')) if '#' in j]
	if len(start_index) != 1:
		print(start_index)
		input("CHECK IT")
	else:
		start_index = start_index[0]

	list_all_data = []

	''' AA, position of AA, chain, secondary structure type, BP1 and BP2, ACC'''
	for stroka in data.strip('\n').split('\n')[start_index+1:]:
		stroka_data = [i for i in stroka.strip('\n').split(' ') if i != '']

		#num_aa,chain,AA,STR
		num_aa = stroka_data[1]
		chain = stroka_data[2]
		aa = stroka_data[3]
		init_secondary_structure = stroka_data[4]

		if aa in SS:
			aa = 'C'

		if not num_aa.isdigit():
			if '-' in num_aa:
				num_aa = num_aa
			elif '!' in num_aa:
				continue
			else:
				chain = num_aa[-1]
				num_aa = num_aa[:-1]
				aa = stroka_data[2]
				init_secondary_structure = stroka_data[3]

		#BP1, BP2, ACC
		bp1, bp2, ACC = find_ACC_bp(stroka_data)

		if bp1 != '0':
			bp1 = '1'

		if bp2 != '0':
			if bp2[0] == '0':
				bp2 = '0'
			else:
				bp2 = '1'

		#Secondary structure type
		if init_secondary_structure in b_sheet:
			secondary_structure = 'B'
		elif init_secondary_structure in a_helix:
			secondary_structure = 'H'
		elif init_secondary_structure in loop:
			secondary_structure = 'O'
		else:
			init_secondary_structure = 'U'
			secondary_structure = 'O'

		list_all_data.append(f'{num_aa},{chain},{aa},{init_secondary_structure},{secondary_structure},{bp1},{bp2},{ACC}')

	''' B-factor '''
	#add b-factor feature
	with open(f"{path_to_read}/{structure}_info.txt", 'r') as file:
		bf_data = file.read().strip('\n').split('\n')

	list_bf_data = []
	filtered_AA = []
	for dssp_stroka in list_all_data:
		check_dssp = dssp_stroka.split(',')[1] + '_' + dssp_stroka.split(',')[2] + '_' + dssp_stroka.split(',')[0]
		for bf_stroka in bf_data:
			check_bf = '_'.join(bf_stroka.split('_')[:-2])
			bf_value = bf_stroka.split('_')[-1]
			if check_dssp == check_bf:
				list_bf_data.append(f"{dssp_stroka},{bf_value}")
				break
		else:
			filtered_AA.append(dssp_stroka)

	if len(filtered_AA) != 0:
		#print(f"{structure_chain}_filter_after_bf.txt is written now!")
		if not os.path.exists(log_path):
			os.mkdir(log_path)
		with open(os.path.join(log_path, "filter_after_bf.txt"), 'w') as file:
			file.write('\n'.join(filtered_AA))

	list_to_dict = lambda x: {'num_aa':[i.split(',')[0] for i in x],
                              'chain':[i.split(',')[1] for i in x],
                              'AA':[i.split(',')[2] for i in x],
                              'init_SS_type':[i.split(',')[3] for i in x],
							  'SS_type':[i.split(',')[4] for i in x],
                              'bp1':[i.split(',')[5] for i in x],
                              'bp2':[i.split(',')[6] for i in x],
                              'ACC':[i.split(',')[7] for i in x],
							  'bfac':[i.split(',')[8] for i in x],
                              }

	result_dict = list_to_dict(list_bf_data)
	result_df = pd.DataFrame(result_dict)

	#change type of some features
	result_df = result_df.astype({'ACC':'int','bfac':'float','bp1':'int','bp2':'int'})

	if len(result_df["bfac"].unique()) == 1:
		if not os.path.exists(log_path):
			os.mkdir(log_path)
		with open(os.path.join(log_path, "constant_bfac.txt"), 'a') as file:
			file.write(f"{structure_chain}\n")

	''' LEN_LOOP '''
	start_index, finish_index = list(result_df.index)[0], list(result_df.index)[-1]

	# Split indices to parts
	def split_indices(indices):
		new_indices = []
		start = 0
		finish = len(indices) - 1
		for i, index in enumerate(indices):
			if i < finish:
				if indices[i+1] - index == 1:
					continue
				else:
					new_indices.append(indices[start:i+1])
					start = i + 1
			else:
				new_indices.append(indices[start:finish+1])
		return new_indices

	indices_B = list(result_df.loc[(result_df["SS_type"] == "B")].index)
	splitted_indices_B = split_indices(indices_B)
	indices_O = list(result_df.loc[(result_df["SS_type"] == "O")].index)
	splitted_indices_O = split_indices(indices_O)
	indices_H = list(result_df.loc[(result_df["SS_type"] == "H")].index)
	splitted_indices_H = split_indices(indices_H)

	result_df["len_loop"] = 0
	for i, j in zip(splitted_indices_O, list(map(len, splitted_indices_O))):
		result_df.loc[i, "len_loop"] = j
	result_df = result_df.astype({"len_loop":"int"})

	""" EDGE_B-STRAND """
	if len(indices_B) == 0:
		result_df["edge_b"] = "0"
	else:
		result_df["temp_feature"] = result_df["bp1"].astype("str") + '_' + result_df["bp2"].astype("str")

		result_df.loc[(result_df["SS_type"] == "B") & (result_df["temp_feature"] == "0_0"), "edge_b"] = "0"
		result_df.loc[(result_df["SS_type"] == "B") & (result_df["temp_feature"] == "1_0"), "edge_b"] = "0"
		result_df.loc[(result_df["SS_type"] == "B") & (result_df["temp_feature"] == "0_1"), "edge_b"] = "0"
		result_df.loc[(result_df["SS_type"] == "B") & (result_df["temp_feature"] == "1_1"), "edge_b"] = "1"

		result_df.loc[(result_df["SS_type"] == "H") | (result_df["SS_type"] == "O"), "edge_b"] = "0"
		result_df = result_df.drop(["temp_feature"], axis=1)

	""" SET THRESHOLDS """
	#ACC_threshold = 200#
	#len_loop_threshold = 20#
	result_df.loc[result_df["ACC"] > 200, "ACC"] = 200
	result_df.loc[result_df["len_loop"] > 20, "len_loop"] = 20

	""" ORDERING AND SAVING """
	if "AF-" in name_file_dssp:
		result_df["structure"] = name_file_dssp.split('.')[0]
	else:
		result_df["structure"] = name_file_dssp.split('_')[0]

	columns = list(result_df.columns)
	result_df = result_df[[columns[-1]] + [columns[1]] + [columns[0]] + columns[2:-1]]

	''' SAVING '''
	result_df.to_csv(f"{path_to_save}/{structure_chain}.csv", index=False)

	''' FINISH '''

""" Set main function """
def main():

	structures = [i for i in os.listdir(script_path) if '.dssp' in i]
	for index, structure in enumerate(structures):
		#print(f"\nStructure: {index+1} --- {structure}")
		form_dataset(script_path, structure, script_path)

		os.remove(os.path.join(script_path, structure))
		os.remove(os.path.join(script_path, f"{structure.split('_')[0]}_seqres.txt"))
		os.remove(os.path.join(script_path, f"{structure.split('_')[0]}_info.txt"))

""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
