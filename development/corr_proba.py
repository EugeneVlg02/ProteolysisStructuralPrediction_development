import os
import time
import glob

import pandas as pd
import plotly.graph_objects as go

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

output_path = os.path.join("/data2/Protease/TEST/output", name_of_set)
if not os.path.exists(output_path):
    os.mkdir(output_path)
save_path = os.path.join(output_path, "corr_proba")
if not os.path.exists(save_path):
        os.mkdir(save_path)

def corr_proba(data, pwm_score_name, structural_score_name, path_to_save):
    name_structure = data["structure"].unique()[0]
    
    colours = {1:{"H":'#006400', "B":'#0000CD', "O":'#FF4500'},
                0:{"H":'#90EE90', "B":'#ADD8E6', "O":'#F0E68C'}}
    symbols = {1:'cross', 0:'circle'}
    names = {1:{"H":{1:'positive helix terminus', 0:'positive helix no terminus'}, "B":{1:'positive b-sheet terminus', 0:'positive b-sheet no terminus'}, "O":{1:'positive loop terminus', 0:'positive loop no terminus'}},
            0:{"H":{1:'negative helix terminus', 0:'negative helix no terminus'}, "B":{1:'negative b-sheet terminus', 0:'negative b-sheet no terminus'}, "O":{1:'negative loop terminus', 0:'negative loop no terminus'}}}
    features = {"ACC_N_chain":True, "bfac_N_chain":False, "len_loop_N_chain":False}

    trace_list = []
    for feature in features.keys():
        for is_cut in [1, 0]:
            for SS_type in ["H", "B", "O"]:
                for is_terminus in [1, 0]:
                    part_of_data = data.loc[(data["is_cut"] == is_cut) & (data["SS_type"] == SS_type) & (data["is_terminus"] == is_terminus)]

                    trace_list.append(go.Scatter(
                        visible=False,
                        x = part_of_data[pwm_score_name],
                        y = part_of_data[structural_score_name],
                        name=names[is_cut][SS_type][is_terminus],
                        mode = "markers",
                        marker = dict(color=colours[is_cut][SS_type], symbol=symbols[is_terminus], size=15 * part_of_data[feature] + 10),
                        customdata = part_of_data[["AA", "num_AA", "SS_type", "init_SS_type", "ACC_N_chain", "bfac_N_chain", "len_loop_N_chain", "is_terminus", "MEROPS_CODE"]],
                        hovertemplate="<b>AA:</b> %{customdata[0]}%{customdata[1]}<br><b>PWM score:</b> %{x:.3f}<br><b>Structural score:</b> %{y:.3f}<br><b>SS_type:</b> %{customdata[3]}:%{customdata[2]}<br><b>ACC (norm):</b> %{customdata[4]:.3f}<br><b>B-factor (norm):</b> %{customdata[5]:.3f}<br><b>Len_loop (norm):</b> %{customdata[6]:.3f}<br><b>Is terminus:</b> %{customdata[7]}<br><b>Protease:</b> %{customdata[8]}"
                        ))

    fig = go.Figure(data=trace_list)
    for i in range(12):
        fig.data[i]["visible"] = True

    fig.update_layout(xaxis={"range":[-35, 15]},
                        xaxis_title="PWM score",
                        yaxis={"range":[0, 1]},
                        yaxis_title="Structural score",
                        title=f"{name_structure}   {pwm_score_name} vs {structural_score_name}",
                        legend=dict(itemsizing="constant")
                        )
    
    # Create slider #
    num_steps = 3
    steps = []
    for i in range(num_steps):
        step = dict(label=list(features.keys())[i],
                    method="restyle",
                    args=["visible", [False]*len(fig.data)])
        for j in range(i * len(fig.data) // 3, (i + 1) * len(fig.data) // 3):
            step["args"][1][j] = True
        steps.append(step)
    slider = [dict(active=0,
                    currentvalue={"prefix":''},
                    steps=steps)]
    fig.layout.sliders = slider

    fig.write_html(path_to_save)

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
            structural_score_column = "LinearDA_Scikit"
            pwm_score_column = [i for i in score_df.columns if "PWM" in i][0]
            save_file = os.path.join(save_path, f"{structure_name}.html")
            corr_proba(score_df, pwm_score_column, structural_score_column, save_file)
            
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")                