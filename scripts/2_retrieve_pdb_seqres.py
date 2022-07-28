""" Import libraries """
import os
import time

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")

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
def get_pdb_seqres(path_to_read, structure_file, path_to_save):
    with open(f"{path_to_read}/{structure_file}", 'r') as file:
        data = file.read()

    seqres_data = [i for i in data.strip('\n').split('\n') if "SEQRES" == i.split(' ')[0]]
    dict_seq = {}

    for seqres in seqres_data:
        clear_seqres = [i for i in seqres.split(' ') if '' != i]

        chain = clear_seqres[2]
        len_chain = clear_seqres[3]
        AA_list = clear_seqres[4:]
        AA_part = ''.join([aminoacid_code[i] if i in aminoacid_code.keys() else "X" for i in AA_list])

        if f"{chain}_{len_chain}" not in dict_seq.keys():
            dict_seq[f"{chain}_{len_chain}"] = AA_part
        else:
            dict_seq[f"{chain}_{len_chain}"] += AA_part

    list_to_write = [f"{i}\n{dict_seq[i]}" for i in dict_seq.keys()]

    with open(f"{path_to_save}/{structure_file.split('.')[0]}_seqres.txt", 'w') as file:
        file.write('\n'.join(list_to_write))

""" Set main function """
def main():
    structures = [i for i in os.listdir(script_path) if ('.pdb' in i) and ('_' not in i)]
    for index, structure in enumerate(structures):
        #print(f"Structure: {index+1} --- {structure}")
        get_pdb_seqres(script_path, structure, script_path)

		#input("Next...")

""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
