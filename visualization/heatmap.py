import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from math import ceil
import sklearn.metrics as metrics
import gzip
from pathlib import Path
import pandas as pd
import json
import seaborn as sns
from matplotlib import rcParams
import gzip
import glob
%matplotlib inline
import csv
# figure size in inches
rcParams['figure.figsize'] = 16,10


root="/projects/p31666/zzhang/doublet-bchmk/final"
files = [f for f in glob.glob(root + "/**/*detection_rates.tsv", recursive=True)]
out_dir="/projects/p31666/zzhang/doublet-bchmk/plots/detection_rate"

doublet_rates_to_test = [0.05, 0.08, 0.1, 0.15, 0.2, 0.25]
for cur_file in files:
    cur_method_id = cur_file.split("/")[6]
    cur_rerun_file = cur_file.replace("final", "partial_rerun")
    cur_res_df = pd.read_csv(cur_file, sep="\t")
    cur_rerun_res_df = pd.read_csv(cur_rerun_file, sep="\t")
    cur_res_all_bchmk_data = cur_res_df.dropna()
    cur_res_df = cur_res_df.dropna()

    # update with rerun data
    cur_res_df.set_index(['ID', 'dbl_exp', 'dbl_act'], inplace=True)
    cur_rerun_res_df.set_index(['ID', 'dbl_exp', 'dbl_act'], inplace=True)
    cur_res_df.update(cur_rerun_res_df)
    cur_res_df.reset_index(inplace=True)
    cur_res_df.to_csv("{}/{}_data.csv".format(out_dir, cur_method_id), index=False)

    # keep as comment
    #     cur_res_df[(cur_res_df["ID"] != "TREX_brain1") &
    #           (cur_res_df["ID"] != "SPLINTR_inVitro_KRAS") &
    #           (cur_res_df["ID"] != "SPLINTR_retransplant")&
    #           (cur_res_df["ID"] != "LARRY_d4_nBC")&
    #           (cur_res_df["ID"] != "LARRY_d4_R_4")]

    all_roc_avg = []
    all_pr_avg = []
    for cur_exp_dbl in doublet_rates_to_test:
        cur_exp_dbl_roc_avg = []
        cur_exp_dbl_pr_avg = []
        for cur_act_dbl in doublet_rates_to_test:
            cur_df = cur_res_all_bchmk_data[(cur_res_all_bchmk_data["dbl_act"]==cur_act_dbl) &
                                            (cur_res_all_bchmk_data["dbl_exp"]==cur_exp_dbl)]
            cur_roc_avg = cur_df["roc"].mean()
            cur_pr_avg = cur_df["pr"].mean()
            cur_exp_dbl_roc_avg.append(cur_roc_avg)
            cur_exp_dbl_pr_avg.append(cur_pr_avg)
        all_roc_avg.append(cur_exp_dbl_roc_avg)
        all_pr_avg.append(cur_exp_dbl_pr_avg)
    plt.clf()
    ax=sns.heatmap(np.array(all_pr_avg), xticklabels=[0.05, 0.08, 0.1, 0.15, 0.2, 0.25], \
                   yticklabels=[0.05, 0.08, 0.1, 0.15, 0.2, 0.25], annot=False, vmin=0.1, vmax=0.7)
    ax.set(xlabel=None, ylabel=None, title=None)
    ax.set(xlabel="", ylabel="")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.xlabel("Actual Doublet Rate")
    plt.ylabel("Expected Doublet Rate")
    plt.title("All Data AUPRC Average: {}".format(cur_method_id))
    plt.savefig("{}/{}_avg_auprc.svg".format(out_dir, cur_method_id))
    plt.clf()

    ax=sns.heatmap(np.array(all_roc_avg), xticklabels=[0.05, 0.08, 0.1, 0.15, 0.2, 0.25], \
                   yticklabels=[0.05, 0.08, 0.1, 0.15, 0.2, 0.25], annot=False, vmin=0.45, vmax=0.9)
    ax.set(xlabel=None, ylabel=None, title=None)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.xlabel("Actual Doublet Rate")
    plt.ylabel("Expected Doublet Rate")
    plt.title("All Data AUROC Average: {}".format(cur_method_id))
    plt.savefig("{}/{}_avg_auroc.svg".format(out_dir, cur_method_id))
    print("{}/{}_avg_auroc.svg".format(out_dir, cur_method_id))
