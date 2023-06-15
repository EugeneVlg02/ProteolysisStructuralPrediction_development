from chimera import runCommand as run
import os
import glob

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = raw_input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

c = 0
for structure in sorted(structure_list):  
    structure_path = os.path.join(StructureSet_path, structure)
    structure_files = glob.glob(os.path.join(structure_path, "*_*.pdb")) + glob.glob(os.path.join(structure_path, "AF-*-F1.pdb"))
    
    for structure_file in structure_files:
        c += 1
        print "{0} --- {1}".format(c, structure_file.split('.')[0].split('/')[-1])
        
        ''' Chimera Run '''
        run("open {}".format(structure_file))
        run("del ligand")
        run("del solvent")
        run("write #0 {}".format(structure_file))
        run("del")

        ''' Processing of output PDB file from chimera: two space must be added! '''
        with open(structure_file, 'r') as file:
            data = file.read()

        new_ATOM = []
        for stroka in data.strip('\n').split('\n'):
            if "ATOM" in stroka or "HETATM" in stroka:
                new_stroka = stroka + ' '*2
                new_ATOM.append(new_stroka)

        with open(structure_file, 'w') as file:
            file.write('\n'.join(new_ATOM))
