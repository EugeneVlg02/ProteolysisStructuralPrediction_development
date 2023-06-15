import os
import time
import glob

import pandas as pd
from chimera import runCommand as run

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = raw_input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

ProCleave_predictions_path = "/data2/Protease/TEST/ProCleave_predictions"

output_path = os.path.join("/data2/Protease/TEST/output", name_of_set)
if not os.path.exists(output_path):
    os.mkdir(output_path)
chimera_path = os.path.join(output_path, "chimera")
if not os.path.exists(chimera_path):
    os.mkdir(chimera_path)
save_scores_path = os.path.join(chimera_path, "map_scores")
if not os.path.exists(save_scores_path):
        os.mkdir(save_scores_path)
save_PWM_path = os.path.join(save_scores_path, "map_PWM2")
if not os.path.exists(save_PWM_path):
        os.mkdir(save_PWM_path)
save_structural_path = os.path.join(save_scores_path, "map_structural2")
if not os.path.exists(save_structural_path):
        os.mkdir(save_structural_path)
save_total_path = os.path.join(save_scores_path, "map_total2")
if not os.path.exists(save_total_path):
        os.mkdir(save_total_path)        
save_ProCleave_path = os.path.join(save_scores_path, "map_ProCleave2")
if not os.path.exists(save_ProCleave_path):
        os.mkdir(save_ProCleave_path)
        
def main():
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        score_path = os.path.join(structure_path, "scores")
        
        score_files = glob.glob(os.path.join(score_path, '*.csv'))
        if len(score_files) == 0: continue
        
        for score_file in score_files:
            structure_name = score_file.split('/')[-1].split('.csv')[0]
            print "{0} --- {1}".format(num, structure_name)
            num += 1
            
            ''' Generate attribute file '''
            score_df = pd.read_csv(score_file, dtype={"structure":str, "chain":str, "num_AA":str, "AA":str, "is_cut":int})
            structure_chain = score_df["structure"].unique()[0]
            
            score_columns = ["LinearDA_Scikit", "Total_score"] + [i for i in score_df.columns if "PWM" in i]
            for score_column in score_columns:
                score_df[score_column + '_attribute'] = '\t:' + score_df["num_AA"] + '.' + score_df["chain"] + '\t' + score_df[score_column].map(str)
                
                if 'Scikit' in score_column:
                    attribute_file = os.path.join(save_structural_path, structure_name + '_' + score_column + '.txt')
                    save_file = os.path.join(save_structural_path, structure_name + '.py')
                elif 'Total' in score_column:
                    attribute_file = os.path.join(save_total_path, structure_name + '_' + score_column + '.txt')
                    save_file = os.path.join(save_total_path, structure_name + '.py')
                else:
                    attribute_file = os.path.join(save_PWM_path, structure_name + '_' + score_column + '.txt')
                    save_file = os.path.join(save_PWM_path, structure_name + '.py')
                    
                with open(attribute_file, 'w') as file:
                    if 'PWM' in score_column:
                        file.write("attribute: PWM\nmatch mode: 1-to-1\nrecipient: residues\n" + '\n'.join(score_df[score_column + "_attribute"].tolist()))
                    else:
                        file.write("attribute: {0}\nmatch mode: 1-to-1\nrecipient: residues\n".format(score_column) + '\n'.join(score_df[score_column + "_attribute"].tolist()))
                ''' Map score '''
                run("del")
                run("open {0}".format(os.path.join(structure_path, structure_chain + '.pdb')))
                run("defattr {}".format(attribute_file))
                
                if 'PWM' in score_column:
                    run("rangecolor PWM -10 blue 2.5 red 15 yellow")
                else:
                    run("rangecolor {0} 0 blue 0.5 red 1 yellow".format(score_column))
                run("save {}".format(save_file))
                os.remove(attribute_file)
            
            ''' ProCleave predictions '''
            ProCleave_df = pd.read_csv(os.path.join(ProCleave_predictions_path, "{0}_ProCleave.csv".format(structure_name)))
            ProCleave_df = ProCleave_df.sort_values("Position").reset_index().rename(columns={"Score":"ProCleave_score", "Position":"num_AA", "Chain Name":"chain"})
            del ProCleave_df["index"]
            ProCleave_df["ProCleave_score_attribute"] = '\t:' + ProCleave_df["num_AA"].map(str) + '.' + ProCleave_df["chain"].map(str) + '\t' + ProCleave_df["ProCleave_score"].map(str)
            
            attribute_file = os.path.join(save_ProCleave_path, structure_name + '_ProCleave.txt')
            save_file = os.path.join(save_ProCleave_path, structure_name + '.py')
            with open(attribute_file, 'w') as file:
                file.write("attribute: ProCleave\nmatch mode: 1-to-1\nrecipient: residues\n" + '\n'.join(ProCleave_df["ProCleave_score_attribute"].tolist()))
            
            run("del")
            run("open {0}".format(os.path.join(structure_path, structure_chain + '.pdb')))
            run("defattr {}".format(attribute_file))
            run("rangecolor ProCleave 0 blue 0.5 red 1 yellow")
            run("save {}".format(save_file))
            os.remove(attribute_file)

if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print("Start: " + str(start) + '\n' + "Finish: " + str(finish) + '\n')