""" Import libraries """
import os
import time
import requests
import json
import argparse

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")
log_path = os.path.join(current_path, "logs")

""" Set URLs """
main_url = 'https://www3.cmbi.umcn.nl/xssp/'

""" Set auxiliary functions """
def pdb_to_dssp(path_to_read, structure_file, main_url):

    url_create = '{}api/create/pdb_file/dssp/'.format(main_url)

    response_create = requests.post(url_create,files={'file_':open(f"{path_to_read}/{structure_file}",'rb')})
    response_create.raise_for_status

    job_id = json.loads(response_create.text)['id']
    #print("Job submitted successfully. Id is: '{}'".format(job_id))

    ready = False
    while not ready:

        url_status = '{}/api/status/pdb_file/dssp/{}/'.format(main_url, job_id)
        response_status = requests.get(url_status)
        response_status.raise_for_status

        status = json.loads(response_status.text)['status']
        #print('Job status is {}'.format(status))

        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE','REVOKED']:
            print('Error')
            break
        else:
            time.sleep(5)
    else:
        url_result = '{}/api/result/pdb_file/dssp/{}/'.format(main_url, job_id)

        response_result = requests.get(url_result)
        response_result.raise_for_status
        result = json.loads(response_result.text)['result']

        print(structure_file.split('.')[0]+'.dssp is created!')
        return result

def pdb_to_dssp_locally(path_to_read,structure_file,path_to_save):
    run_dssp = f"mkdssp -i {path_to_read}/{structure_file} -o {path_to_save}/{structure_file.split('.')[0]}.dssp"
    os.system(run_dssp)

""" Set main function """
def main():
    parser = argparse.ArgumentParser(description='### Here should be the description! ###')
    parser.add_argument("-l", "--launch", help="Launch DSSP locally (specify '-l') or via server", action="store_true")
    args = parser.parse_args()

    error_dssp = []

    structures = [i for i in os.listdir(script_path) if ('.pdb' in i) and ('_' in i)]
    for index, structure in enumerate(structures):
        #print(f"\nStructure: {index+1} --- {structure}")

        if args.launch:
            pdb_to_dssp_locally(script_path, structure, script_path)
            if not os.path.exists(os.path.join(script_path, f"{structure.split('.')[0]}.dssp")):
                error_dssp.append(f"{structure} - error DSSP! (locally)")
        else:
            result = pdb_to_dssp(script_path, structure, main_url)
            if result != None:
                with open(os.path.join(script_path, f"{structure.split('.')[0]}.dssp"), 'w') as file:
                    file.write(result)
            else:
                error_dssp.append(f"{structure} --- error DSSP! (via server)")

    if len(error_dssp) > 0:
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        with open(os.path.join(log_path, "not_dssp_files.txt"), 'w') as file:
            file.write('\n'.join(error_dssp))

    os.remove(os.path.join(script_path, structure))
    os.remove(os.path.join(script_path, f"{structure.split('_')[0]}.pdb"))

""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
