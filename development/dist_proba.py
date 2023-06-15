import os
import glob
import time

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

output_path = os.path.join("/data2/Protease/TEST/output", name_of_set)
if not os.path.exists(output_path):
    os.mkdir(output_path)
save_path = os.path.join(output_path, "dist_proba")
if not os.path.exists(save_path):
        os.mkdir(save_path)

def main():
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        score_path = os.path.join(structure_path, "scores")
        
        score_files = glob.glob(os.path.join(score_path, '*.csv'))
        if len(score_files) == 0: continue
        
        for score_file in score_files:
            structure_name = score_file.split('/')[-1].split('.csv')[0]
            print(f"{num} --- {structure_name}")
            num += 1
            
            score_df = pd.read_csv(score_file, dtype={"structure":str, "chain":str, "num_AA":str, "AA":str})
            
            scores = score_df["LinearDA_Scikit"].tolist() + score_df["Total_score"].tolist()
            labels = ["Structural score"] * len(score_df) + ["Total score"] * len(score_df)
            plot_df = pd.DataFrame({"Score":scores, "Type":labels})
            
            fig = make_subplots(rows=1, cols=3, subplot_titles=("Structural score", "Total score", "PWM score"))
            fig.add_trace(
                          go.Box(y=score_df["LinearDA_Scikit"].tolist(), name="Structural score"),
                          row=1, col=1
                          )
            fig.add_trace(
                          go.Box(y=score_df["Total_score"].tolist(), name="Total score"),
                          row=1, col=2
                          )
            fig.add_trace(
                          go.Box(y=score_df[[i for i in score_df.columns if "PWM" in i][0]].tolist(), name="PWM score"),
                          row=1, col=3
                          )
            
            fig.update_xaxes(showticklabels=False)
            fig.update_yaxes(title_text="Score", range=[0, 1], row=1, col=1)
            fig.update_yaxes(range=[0, 1], row=1, col=2)
            
            fig.update_layout(title_text=structure_name)
            fig.write_html(os.path.join(save_path, f"{structure_name}.html"))

''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")