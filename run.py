
import shutil
import os
import time
import argparse

start = time.ctime()

current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")

def main():
    print("######## START ########")

    parser = argparse.ArgumentParser(description='### Here should be the description! ###')
    parser.add_argument("-input", help="The name of your PDB ID with separate chain, for example '4GAW_A'", type=str)
    parser.add_argument("-model", help="The name of model used, for example 'GaussianNB'", type=str)
    parser.add_argument("-protease", help="The name of protease used, for example 'MMP9'", type=str)
    parser.add_argument("-pwm", help="The name of PWM used, for example 'PMAP'", type=str)
    parser.add_argument("-l", "--launch", help="Launch DSSP locally (specify '-l') or via server", action="store_true")
    args = parser.parse_args()

    pdb_id = args.input
    model = args.model
    protease = args.protease
    pwm = args.pwm

    print("######## Step 1 -- Creating of PDB file with separate chain ########")
    os.system(f"python {script_path}/1_download_structure_with_sep_chain.py -i {pdb_id}")

    print("######## Step 2 -- Preprocessing  ########")
    os.system(f"python {script_path}/2_retrieve_pdb_seqres.py")
    os.system(f"python {script_path}/3_retrieve_pdb_atom.py")

    print("######## Step 3 -- Creating the DSSP file ########")
    if args.launch:
        os.system(f"python {script_path}/4_create_dssp.py -l")
    else:
        os.system(f"python {script_path}/4_create_dssp.py")

    print("######## Step 4 -- Creating the ad-hoc features ########")
    os.system(f"python {script_path}/5_create_adhoc.py")
    os.system(f"python {script_path}/6_add_edge_parts.py")

    print("######## Step 5 -- Normalisation of data ########")
    os.system(f"python {script_path}/7_norm_data.py")

    print("######## Step 6 -- Predicting and visualising of scores ########")
    os.system(f"python {script_path}/8_get_scores.py -input {pdb_id} -model {model} -protease {protease} -pwm {pwm}")
    os.system(f"python {script_path}/9_visualise_scores.py -input {pdb_id} -model {model} -protease {protease} -pwm {pwm}")

    print("######## FINISH ########")
if __name__ == "__main__":

    main()

    end = time.ctime()
    print(f"\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n")
