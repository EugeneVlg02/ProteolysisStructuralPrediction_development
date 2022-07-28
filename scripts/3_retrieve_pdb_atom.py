""" Import libraries """
import os
import time
import functools

import numpy as np

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")
log_path = os.path.join(current_path, "logs")

""" Set AA code """
aminoacid_code = {'GLY':'G',
				  'ALA':'A',
				  'VAL':'V',
				  'LEU':'L',
				  'ILE':'I',
				  'CYS':'C',
				  'MET':'M',
				  'PHE':'F',
				  'TYR':'Y',
				  'TRP':'W',
				  'PRO':'P',
				  'SER':'S',
				  'THR':'T',
				  'ASN':'N',
				  'GLN':'Q',
				  'ASP':'D',
				  'GLU':'E',
				  'HIS':'H',
				  'LYS':'K',
				  'ARG':'R'}

""" Set auxiliary function """
def retrieve_data_from_atom(atom_stroka):

    if not atom_stroka[0].isalpha():
        possible_AA = atom_stroka[2]
        possible_chain = atom_stroka[3]
        possible_position = atom_stroka[4]
    elif len(atom_stroka) == 13:
        possible_AA = atom_stroka[4]
        possible_chain = atom_stroka[5]
        possible_position = atom_stroka[6]
    else:
        possible_AA = atom_stroka[3]
        possible_chain = atom_stroka[4]
        possible_position = atom_stroka[5]

    bfac = atom_stroka[-2]

    if '.' in possible_chain:
        AA = atom_stroka[2][-3:]
        chain = possible_AA[0]
        if '-' in possible_AA[1:]:
            chain = f"{chain}+-"
            position = possible_AA[2:]
        else:
            if not possible_AA[1:].isdigit():
                chain = f"{chain}+{possible_AA[-1]}"
                position = possible_AA[1:-1]
            else:
                position = possible_AA[1:]
    else:
        if '.' in possible_position:
            if atom_stroka[2][-3:] in aminoacid_code.keys():
                AA = atom_stroka[2][-3:]
                chain = possible_AA
                if '-' in possible_chain:
                    chain = f"{chain}+-"
                    position = possible_chain[1:]
                else:
                    if not possible_chain.isdigit():
                        chain = f"{chain}+{possible_chain[-1]}"
                        position = possible_chain[:-1]
                    else:
                        position = possible_chain
            else:
                AA = possible_AA[-3:]
                chain = possible_chain[0]
                if '-' in possible_chain[1:]:
                    chain = f"{chain}+-"
                    position = possible_chain[2:]
                else:
                    if not possible_chain[1:].isdigit():
                        chain = f"{chain}+{possible_chain[-1]}"
                        position = possible_chain[1:-1]
                    else:
                        position = possible_chain[1:]
        else:
            AA = possible_AA[-3:]
            chain = possible_chain
            if '-' in possible_position:
                chain = f"{chain}+-"
                position = possible_position[1:]
            else:
                if not possible_position.isdigit():
                    chain = f"{chain}+{possible_position[-1]}"
                    position = possible_position[:-1]
                else:
                    position = possible_position

        # 1WAR A 76 #
    if bfac.count(".") == 2:
        bfac = bfac[4:]

    return AA, chain, position, bfac

""" Set function for analysing of structure PDB file """
def pdb_file_analysis(path_to_read, structure_file, path_to_save):
	with open(f"{path_to_read}/{structure_file}", 'r') as file:
		data = file.read()

	TER_string = [(i,j) for i,j in enumerate(data.strip('\n').split('\n')) if j.split(' ')[0] == 'TER']
	end_TER_index = TER_string[-1][0]

	missing_res = [i for i in data.strip('\n').split('\n') if (i.split(' ')[0] == "REMARK") and (i.split(' ')[1]) == "465"]
	start_index_mis = [i for i,j in enumerate(missing_res) if "SSSEQI" in j]

	if len(start_index_mis) == 1:
		start_index_mis = start_index_mis[0]
	elif len(start_index_mis) == 0:
		start_index_mis = -1
	else:
		print(start_index_mis)
		input("CHECK IT")

	clear_missing_res = [[j for j in i.split(' ') if '' != j] for i in missing_res][start_index_mis+1:]

	atom = [i for j,i in enumerate(data.strip('\n').split('\n')) if (i.split(' ')[0] == "ATOM") or ("HETATM" in i.split(' ')[0] and j < end_TER_index)]
	clear_atom = [[j for j in i.split(' ') if '' != j] for i in atom]

	dict_atom = {}
	dict_bf = {}
	for stroka in clear_atom:
		atom_AA, atom_chain, atom_position, atom_bfac = retrieve_data_from_atom(stroka)

		if atom_AA in aminoacid_code.keys():
			if '+' in atom_chain:
				if '-' in atom_chain:
					dict_atom.setdefault(atom_chain.split('+')[0], [])
					dict_bf.setdefault(f"{aminoacid_code[atom_AA]}_{atom_chain.split('+')[0]}_-{atom_position}", [])
					if f"{aminoacid_code[atom_AA]}_-{atom_position}" not in dict_atom[atom_chain.split('+')[0]]:
						dict_atom[atom_chain.split('+')[0]].append(f"{aminoacid_code[atom_AA]}_-{atom_position}")
					if float(atom_bfac) not in dict_bf[f"{aminoacid_code[atom_AA]}_{atom_chain.split('+')[0]}_-{atom_position}"]:
						dict_bf[f"{aminoacid_code[atom_AA]}_{atom_chain.split('+')[0]}_-{atom_position}"].append(float(atom_bfac))
				else:
					dict_atom.setdefault(atom_chain.split('+')[0], [])
					dict_bf.setdefault(f"{aminoacid_code[atom_AA]}_{atom_chain.split('+')[0]}_{atom_position}{atom_chain.split('+')[1]}", [])
					if f"{aminoacid_code[atom_AA]}_{atom_position}{atom_chain.split('+')[1]}" not in dict_atom[atom_chain.split('+')[0]]:
						dict_atom[atom_chain.split('+')[0]].append(f"{aminoacid_code[atom_AA]}_{atom_position}{atom_chain.split('+')[1]}")
					if float(atom_bfac) not in dict_bf[f"{aminoacid_code[atom_AA]}_{atom_chain.split('+')[0]}_{atom_position}{atom_chain.split('+')[1]}"]:
						dict_bf[f"{aminoacid_code[atom_AA]}_{atom_chain.split('+')[0]}_{atom_position}{atom_chain.split('+')[1]}"].append(float(atom_bfac))
			else:
				dict_atom.setdefault(atom_chain, [])
				dict_bf.setdefault(f"{aminoacid_code[atom_AA]}_{atom_chain}_{atom_position}", [])
				if f"{aminoacid_code[atom_AA]}_{atom_position}" not in dict_atom[atom_chain]:
					dict_atom[atom_chain].append(f"{aminoacid_code[atom_AA]}_{atom_position}")
				if float(atom_bfac) not in dict_bf[f"{aminoacid_code[atom_AA]}_{atom_chain}_{atom_position}"]:
					dict_bf[f"{aminoacid_code[atom_AA]}_{atom_chain}_{atom_position}"].append(float(atom_bfac))
		else:
			if '+' in atom_chain:
				if '-' in atom_chain:
					dict_atom.setdefault(atom_chain.split('+')[0], [])
					dict_bf.setdefault(f"X_{atom_chain.split('+')[0]}_-{atom_position}", [])
					if f"X_-{atom_position}" not in dict_atom[atom_chain.split('+')[0]]:
						dict_atom[atom_chain.split('+')[0]].append(f"X_-{atom_position}")
					if float(atom_bfac) not in dict_bf[f"X_{atom_chain.split('+')[0]}_-{atom_position}"]:
						dict_bf[f"X_{atom_chain.split('+')[0]}_-{atom_position}"].append(float(atom_bfac))
				else:
					dict_atom.setdefault(atom_chain.split('+')[0], [])
					dict_bf.setdefault(f"X_{atom_chain.split('+')[0]}_{atom_position}{atom_chain.split('+')[1]}", [])
					if f"X_{atom_position}{atom_chain.split('+')[1]}" not in dict_atom[atom_chain.split('+')[0]]:
						dict_atom[atom_chain.split('+')[0]].append(f"X_{atom_position}{atom_chain.split('+')[1]}")
					if float(atom_bfac) not in dict_bf[f"X_{atom_chain.split('+')[0]}_{atom_position}{atom_chain.split('+')[1]}"]:
						dict_bf[f"X_{atom_chain.split('+')[0]}_{atom_position}{atom_chain.split('+')[1]}"].append(float(atom_bfac))
			else:
				dict_atom.setdefault(atom_chain.split('+')[0], [])
				dict_bf.setdefault(f"X_{atom_chain}_{atom_position}", [])
				if f"X_{atom_position}" not in dict_atom[atom_chain]:
					dict_atom[atom_chain].append(f"X_{atom_position}")
				if float(atom_bfac) not in dict_bf[f"X_{atom_chain}_{atom_position}"]:
					dict_bf[f"X_{atom_chain}_{atom_position}"].append(float(atom_bfac))

	dict_atom = {i:[f"{j}_str" for j in dict_atom[i]] for i in dict_atom.keys()}
	dict_bf = {i:str(round(np.mean(dict_bf[i]),2)) for i in dict_bf.keys()}

	dict_missing_residues = {}

	for stroka in clear_missing_res:
		if len(stroka) == 5:
			mis_AA, mis_chain, mis_position = stroka[2], stroka[3], stroka[4]
			mis_model = "1"
		else:
			mis_model, mis_AA, mis_chain, mis_position = stroka[2], stroka[3], stroka[4], stroka[5]

		if mis_model != "1": continue

		dict_missing_residues.setdefault(mis_chain, [])

		if mis_AA in aminoacid_code.keys():
			if f"{aminoacid_code[mis_AA]}_{mis_position}" not in dict_missing_residues[mis_chain]:
				dict_missing_residues[mis_chain].append(f"{aminoacid_code[mis_AA]}_{mis_position}")
		else:
			if f"X_{mis_position}" not in dict_missing_residues[mis_chain]:
				dict_missing_residues[mis_chain].append(f"X_{mis_position}")

	dict_missing_residues = {i:[f"{j}_mis" for j in dict_missing_residues[i]] for i in dict_missing_residues.keys()}

	def my_sort_function(a,b):
		if a.split('_')[1][-1].isalpha():
			if b.split('_')[1][-1].isalpha():
				if int(a.split('_')[1][:-1]) > int(b.split('_')[1][:-1]):
					return 1
				elif int(a.split('_')[1][:-1]) < int(b.split('_')[1][:-1]):
					return -1
				else:
					if a.split('_')[1][-1] > b.split('_')[1][-1]:
						return 1
					else:
						return -1
			else:
				if int(a.split('_')[1][:-1]) >= int(b.split('_')[1]):
					return 1
				else:
					return -1
		else:
			if b.split('_')[1][-1].isalpha():
				if int(a.split('_')[1]) > int(b.split('_')[1][:-1]):
					return 1
				else:
					return -1
			else:
				if int(a.split('_')[1]) > int(b.split('_')[1]):
					return 1
				else:
					return -1

	#SEQRES#
	with open(f"{path_to_read}/{structure_file.split('.')[0]}_seqres.txt", 'r') as file:
		data = file.read()
	dict_seqres = {j:data.strip('\n').split('\n')[i+1] for i,j in enumerate(data.strip('\n').split('\n')) if i%2 == 0}
	chains = [i.split('_')[0] for i in dict_seqres.keys()]

	common_dict = {}
	for chain in chains:
		if chain not in dict_missing_residues.keys():
			common_dict[chain] = dict_atom[chain]
		else:
			common_dict[chain] = dict_atom[chain] + dict_missing_residues[chain]
			common_dict[chain].sort(key=functools.cmp_to_key(my_sort_function))

	list_to_write = []
	for chain_len in dict_seqres.keys():
		chain,len_chain = chain_len.split('_')[0],chain_len.split('_')[1]

		ATOM_seq = ''.join([i.split('_')[0] for i in common_dict[chain]])
		SEQRES_seq = dict_seqres[chain_len]
		if ATOM_seq == SEQRES_seq:
			list_to_write.append('\n'.join([f"{chain}_{i}" for i in common_dict[chain]]))
		else:
			#print(chain)#
			#print(ATOM_seq)#
			#print(SEQRES_seq)#
			#input("NOT CORRESPOND")#
			if not os.path.exists(log_path):
				os.mkdir(log_path)
			with open(os.path.join(log_path, "not_correspond_index_aa.txt"), 'a') as file:
				file.write(structure_file.split('.')[0] + '\t' + chain + '\n')

	list_to_write_2 = ['\n'.join([j+'_'+dict_bf[j.split('_')[1]+'_'+j.split('_')[0]+'_'+j.split('_')[2]] if j.split('_')[1]+'_'+j.split('_')[0]+'_'+j.split('_')[2] in dict_bf.keys() else j for j in i.split('\n')]) for i in list_to_write]

	with open(f"{path_to_read}/{structure_file.split('.')[0]}_info.txt", 'w') as file:
		file.write('\n'.join(list_to_write_2))

""" Set main function """
def main():
	structures = [i for i in os.listdir(script_path) if ('.pdb' in i) and ('_' not in i)]
	for index, structure in enumerate(structures):
		#print(f"Structure: {index+1} --- {structure}")
		pdb_file_analysis(script_path, structure, script_path)

""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
