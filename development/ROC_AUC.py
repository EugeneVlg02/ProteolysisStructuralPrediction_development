''' Import libraries '''
import os
import glob
import time

from sklearn.metrics import roc_auc_score, roc_curve
import pandas as pd
import numpy as np
import plotly.express as px

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

ProCleave_predictions_path = "/data2/Protease/TEST/ProCleave_predictions"
output_path = os.path.join("/data2/Protease/TEST/output", name_of_set)
if not os.path.exists(output_path):
    os.mkdir(output_path)
save_path = os.path.join(output_path, "ROC_curves")
if not os.path.exists(save_path):
        os.mkdir(save_path)
output_path = os.path.join("/data2/Protease/TEST/output", "TEST_MEROPS")

def main():
    total_mean = []
    str_mean = []
    PWM_mean = []
    ProCleave_mean = []
    total_cuts = 0
    
    test_set = pd.DataFrame()
    compare_set = pd.DataFrame()
    
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
            if structure_name in ["1AUI_A_C14.003", "5KKR_C_M10.003"]:
                print("MK exclude")
                continue
            
            score_df = pd.read_csv(score_file, dtype={"structure":str, "chain":str, "num_AA":str, "AA":str})
            test_set = pd.concat([test_set, score_df], ignore_index=True)
            
            score_df = score_df.iloc[3:-3].reset_index()
            del score_df["index"]
            ProCleave_df = pd.read_csv(os.path.join(ProCleave_predictions_path, f"{structure_name}_ProCleave.csv"))
            ProCleave_df = ProCleave_df.sort_values("Position").reset_index().rename(columns={"Score":"ProCleave_score"})
            del ProCleave_df["index"]
            
            compare_df = pd.concat([score_df, ProCleave_df["ProCleave_score"]], axis=1)
            compare_set = pd.concat([compare_set, compare_df], ignore_index=True)
            
            num_cuts = compare_df.loc[compare_df["is_cut"] == 1, "is_cut"].sum()
            total_cuts += num_cuts
            
            ''' ROC AUC and ROC-curves '''
            
            true = compare_df["is_cut"].tolist()
            total_proba = compare_df["GaussianNaiveBayes_Total_score"].tolist()
            str_proba = compare_df["LinearDA_Scikit"].tolist()
            PWM_proba = compare_df[[i for i in compare_df.columns if 'PWM' in i][0]].tolist()
            ProCleave_proba = compare_df["ProCleave_score"].tolist()
            
            total_curve = roc_curve(true, total_proba)
            total_AUC = round(roc_auc_score(true, total_proba), 3)
            total_mean.append(total_AUC)
            total_df = pd.DataFrame({"FPR":total_curve[0], "TPR":total_curve[1], "Type":f"Total score with AUC: {total_AUC}"})
            
            str_curve = roc_curve(true, str_proba)
            str_AUC = round(roc_auc_score(true, str_proba), 3)
            str_mean.append(str_AUC)
            str_df = pd.DataFrame({"FPR":str_curve[0], "TPR":str_curve[1], "Type":f"Structural score with AUC: {str_AUC}"})
            
            PWM_curve = roc_curve(true, PWM_proba)
            PWM_AUC = round(roc_auc_score(true, PWM_proba), 3)
            PWM_mean.append(PWM_AUC)
            PWM_df = pd.DataFrame({"FPR":PWM_curve[0], "TPR":PWM_curve[1], "Type":f"PWM score with AUC: {PWM_AUC}"})
            
            ProCleave_curve = roc_curve(true, ProCleave_proba)
            ProCleave_AUC = round(roc_auc_score(true, ProCleave_proba), 3)
            ProCleave_mean.append(ProCleave_AUC)
            ProCleavePred_df = pd.DataFrame({"FPR":ProCleave_curve[0], "TPR":ProCleave_curve[1], "Type":f"ProCleave score with AUC: {ProCleave_AUC}"})
                  
            plot_df = pd.concat([total_df, str_df, PWM_df, ProCleavePred_df], ignore_index=True)
            fig = px.line(plot_df, x="FPR", y="TPR", color="Type", title=f"{structure_name}    The number of proteolytic events: {num_cuts}")
            #fig.write_html(os.path.join(save_path, f"{structure_name}_ProCleave.html"))
    
    plot_total_df = pd.DataFrame({"Method":"CutProba", "AUC":total_mean})
    plot_ProCleave_df = pd.DataFrame({"Method":"ProCleave", "AUC":ProCleave_mean})
    plot_df = pd.concat([plot_total_df, plot_ProCleave_df], ignore_index=True)
    print(plot_df)
    #plot_df.to_csv(os.path.join(output_path, f"CutProba_vs_ProCleave_TEST_MEROPS.csv"), index=False)
    
    fig = px.box(plot_df, x="Method", y="AUC", color="Method", width=650)
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(range=[0.45, 1.15], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>Test set</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="Method",
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                boxgap=0.5,
                )
    #fig.write_html(os.path.join(output_path, f"CutProba_vs_ProCleave_TEST_MEROPS.html"))
    #fig.write_image(os.path.join(output_path, f"CutProba_vs_ProCleave_TEST_MEROPS.jpg"), width = 650, height = 550, scale = 4)
    print(np.median(total_mean))
    print(np.median(ProCleave_mean))
    print(f"The number of test structures: {len(total_mean)}")
    print(f"The total number of cuts: {total_cuts}")
    print(f"The mean of total AUC: {round(np.mean(total_mean), 3)}")
    print(f"The mean of structural AUC: {round(np.mean(str_mean), 3)}")
    print(f"The mean of PWM AUC: {round(np.mean(PWM_mean), 3)}")
    print(f"The mean of ProCleave AUC: {round(np.mean(ProCleave_mean), 3)}")
    test_set.to_csv(os.path.join(output_path, "TEST_MEROPS.csv"), index=False)
    compare_set.to_csv(os.path.join(output_path, "Our_vs_ProCleave.csv"), index=False)
    
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
