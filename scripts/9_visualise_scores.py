""" Import libraries """
import os
import time
import argparse

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")

import plotly
import plotly.graph_objs as go

""" Set paths """
current_path = os.getcwd()
fitted_models_path = os.path.join(current_path, "fitted_models")
pwmData_path = os.path.join(current_path, "PWM")
script_path = os.path.join(current_path, "scripts")

parser = argparse.ArgumentParser(description='### Here should be the description! ###')
parser.add_argument("-input", help="The name of your PDB ID with separate chain, for example '4GAW_A'", type=str)
parser.add_argument("-model", help="The name of model used, for example 'GaussianNB'", type=str)
parser.add_argument("-protease", help="The name of protease used, for example 'MMP9'", type=str)
parser.add_argument("-pwm", help="The name of PWM used, for example 'PMAP'", type=str)

args = parser.parse_args()

structure_chain = args.input
model_name = args.model
protease_name = args.protease
pwm_name = args.pwm

feature_path = os.path.join(current_path, structure_chain)

pwm_score_name = f"{protease_name}_{pwm_name}_score"
model_score_name = f"{model_name}_score"

col_type = {"structure":str, "chain":str, "num_aa":str, "substrate":str}

# Auxiliary function 2 #
def visualise_sns(data, pwm_score_name, structural_score_name, path_to_save):
    fig, ax = plt.subplots(1, 1, figsize=(19, 10))
    sns.scatterplot(x=pwm_score_name, y=structural_score_name, data=data, s=100, ax=ax)
    ax.set_xlabel("PWM score", fontsize=14)
    ax.set_ylabel("Structural score", fontsize=14)
    ax.set_xlim([-35, 15])
    ax.set_ylim([0, 1])
    ax.set_title(f"{pwm_score_name.split('_score')[0]} vs {structural_score_name.split('_')[0]}", fontsize=18, fontweight="bold")
    plt.tight_layout()

    plt.savefig(path_to_save)

def visualise_plotly(data, pwm_score_name, structural_score_name, path_to_save):
    colours = {"H":'#006400', "B":'#0000CD', "O":'#FF4500'}
                #0:{"H":'#90EE90', "B":'#ADD8E6', "O":'#F0E68C'}}
    symbols = {1:'cross', 0:'circle'}
    names = {"H":{1:'helix edge', 0:'helix noedge'}, "B":{1:'b-sheet edge', 0:'b-sheet noedge'}, "O":{1:'loop edge', 0:'loop noedge'}}
            #{1:{"H":{1:'positive helix edge', 0:'positive helix noedge'}, "B":{1:'positive b-sheet edge', 0:'positive b-sheet noedge'}, "O":{1:'positive loop edge', 0:'positive loop noedge'}},
    features = {"ACC_N_com":True, "bfac_N_chain":False, "len_loop_N_chain":False}

    trace_list = []
    for feature in features.keys():
        for SS_type in ["H", "B", "O"]:
            for is_edge in [1, 0]:
                part_of_data = data.loc[(data["SS_type"] == SS_type) & (data["is_edge"] == is_edge)]

                trace_list.append(go.Scatter(
                    visible=False,
                    x = part_of_data[pwm_score_name],
                    y = part_of_data[structural_score_name],
                    name=names[SS_type][is_edge],
                    mode = "markers",
                    marker = dict(color=colours[SS_type], symbol=symbols[is_edge], size=15*part_of_data[feature]+10),
                    customdata = part_of_data[["AA", "num_aa", "SS_type", "init_SS_type", "ACC_N_com", "bfac_N_chain", "len_loop_N_chain", "is_edge"]],
                    hovertemplate="<b>AA:</b> %{customdata[0]}%{customdata[1]}<br><b>PWM score:</b> %{x:.3f}<br><b>Structural score:</b> %{y:.3f}<br><b>SS_type:</b> %{customdata[3]}:%{customdata[2]}<br><b>ACC (norm):</b> %{customdata[4]:.3f}<br><b>B-factor (norm):</b> %{customdata[5]:.3f}<br><b>Len_loop (norm):</b> %{customdata[6]:.3f}<br><b>Is edge:</b> %{customdata[7]}"
                    ))

    fig = go.Figure(data=trace_list)
    for i in range(6):
        fig.data[i]["visible"] = True

    fig.update_layout(xaxis={"range":[-35, 15]},
                        xaxis_title="PWM score",
                        yaxis={"range":[0, 1]},
                        yaxis_title="Structural score",
                        title=f"{pwm_score_name} vs {structural_score_name}",
                        legend=dict(itemsizing="constant")
                        )
    # Create slider #

    num_steps = 3
    steps = []
    for i in range(num_steps):
        step = dict(label=list(features.keys())[i],
                    method="restyle",
                    args=["visible", [False]*len(fig.data)])
        for j in range(i*len(fig.data)//3, (i+1)*len(fig.data)//3):
            step["args"][1][j] = True
        steps.append(step)
    slider = [dict(active=0,
                    currentvalue={"prefix":''},
                    steps=steps)]
    fig.layout.sliders = slider

    fig.write_html(path_to_save)

def main():

    score_files = [i for i in os.listdir(feature_path) if ('.csv' in i) and (protease_name in i) and (pwm_name in i) and (model_name in i)]
    for score_file in score_files:
        score_df = pd.read_csv(os.path.join(feature_path, score_file), dtype=col_type)

        savefile_name1 = score_file.split('.')[0] + ".html"
        savefile_name2 = score_file.split('.')[0] + ".jpg"
        visualise_plotly(score_df, pwm_score_name, model_score_name, os.path.join(feature_path, savefile_name1))
        visualise_sns(score_df, pwm_score_name, model_score_name, os.path.join(feature_path, savefile_name2))


""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
