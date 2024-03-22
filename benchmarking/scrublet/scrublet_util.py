import scrublet as scr
import scipy.io
import numpy as np
import os
import sklearn.metrics as metrics
import gzip
import pandas as pd
from scipy.sparse import vstack
import random
import re


def run_scrublet_bchmk_data(data_dir):
    # ValueError: n_components=30 must be between 1 and min(n_samples, n_features)=9 with svd_solver='arpack'
    cur_scrub, doublet_labels = load_bchmk_into_scrublet(data_dir)
    doublet_scores, predicted_doublets = cur_scrub.scrub_doublets()
    fpr, tpr, threshold = metrics.roc_curve(y_true=doublet_labels, y_score=doublet_scores, pos_label="doublet")
    roc_auc = metrics.auc(fpr, tpr)
    precision, recall, thresholds = metrics.precision_recall_curve(doublet_labels, doublet_scores, pos_label="doublet")
    auc_precision_recall = metrics.auc(recall, precision)
    out_df = pd.DataFrame({
        "roc": [roc_auc],
        "pr": [auc_precision_recall]
    })
    df = pd.DataFrame({'score': doublet_scores, 'label': doublet_labels, "predictions": predicted_doublets})
    return out_df, df


def run_bchmk_on_bchmk_data(data_dir, data_sample_ID, num_to_test):
    # save_dir = "/projects/b1042/GoyalLab/zzhang/doublet_objects_benchmarking_data/scrublet"
    save_dir = "/projects/b1042/GoyalLab/zzhang/doublet_objects_benchmarking_data_no_mod/scrublet"
    df_ls = []
    for idx in range(num_to_test):
        cur_res, df = run_scrublet_bchmk_data(data_dir)
        df_ls.append(cur_res)
        outfile = f"{save_dir}/{data_sample_ID}.csv"
        df.to_csv(outfile, index = False)

    all_stats = pd.concat(df_ls)
    all_stats['ID'] = data_sample_ID
    return all_stats


def load_bchmk_into_scrublet(data_dir):
    dataset_id = data_dir.split("/")[-1]
    # count_file = data_dir + "/" + dataset_id + "_0.08.mtx"
    # labels_file = data_dir + "/" + dataset_id + "_labels_0.08.csv"
    count_file = data_dir + "/" + dataset_id + ".mtx"
    print(count_file)
    labels_file = data_dir + "/" + dataset_id + "_labels.csv"
    labels = list(pd.read_csv(labels_file, header=None)[0])
    count_mtx = scipy.io.mmread(count_file).T
    cur_scrub = scr.Scrublet(counts_matrix=count_mtx)
    return cur_scrub, labels


def load_fatemap_into_scrublet(data_dir, exp_dbl, doublets_pct=0.08, dbl_method = "sum"):
    cell_id_file = data_dir + "/barcodes.tsv.gz"
    gene_id_file = data_dir + "/features.tsv.gz"
    mtx = data_dir + "/matrix.mtx.gz"
    singlets_file = data_dir + "/singlets_all.txt"
    singlets_df = pd.read_csv(singlets_file, header=None)
    singlets_pair_for_doublet_simulation_file = data_dir + "/corrected_singlet_pairs.csv"
    singlets_pairs_df = pd.read_csv(singlets_pair_for_doublet_simulation_file, header=None, skipinitialspace=True)
    singlets_pairs_df[0] = singlets_pairs_df[0]

    with gzip.open(cell_id_file) as fp:
        cell_ids = np.array([x.decode('UTF-8').strip("\n") for x in fp.readlines()])

    # check if "-1" is a part of 10X barcodes
    if not ("watermelon" in data_dir):
        if "-1" in cell_ids[0]:
            singlets_df[0] = singlets_df[0] + "-1"
        else:
            singlets_df[0] = singlets_df[0].str.replace("-1", "", regex=False)
    else:
        plus_one_str = re.search(r'-(\d+)$', cell_ids[0]).group(0)
        singlets_df[0] = singlets_df[0] + plus_one_str

    # identify singlets present in 10X
    singlets_in_seurat = [x for x in singlets_df[0] if x in cell_ids]
    # get index of singlets_in_seurat, in its order
    # basically now the count matrix is in the order of singlets in seurat
    idx_to_keep = [np.where(cell_ids == x)[0][0] for x in singlets_in_seurat]
    # cell id by gene id matrix
    # in the order of singlets in seurat
    counts_matrix = scipy.io.mmread(mtx).T.tocsr()[idx_to_keep, :]
    #     print("Total singlets: {}".format(str(counts_matrix.shape[0])))

    # take average of expression to produce doublets
    singlets_pt1_idx = [singlets_in_seurat.index(x) for x in singlets_pairs_df[0]]
    singlets_pt2_idx = [singlets_in_seurat.index(x) for x in singlets_pairs_df[1]]
    if dbl_method == "avg":
        print("INFO: avg doublets!")
        avg_assay_ori = (counts_matrix[singlets_pt1_idx, :] + counts_matrix[singlets_pt2_idx, :]) / 2
    else:
        print("INFO: summation doublets!")
        avg_assay_ori = (counts_matrix[singlets_pt1_idx, :] + counts_matrix[singlets_pt2_idx, :])
    doublet_id = pd.DataFrame(singlets_pairs_df[0] + "--" + singlets_pairs_df[1])
    doublet_id = doublet_id[0].str.strip()
    num_doublets = int((len(singlets_in_seurat) * doublets_pct) / (1 + doublets_pct))
    avg_assay = avg_assay_ori[:num_doublets, :]

    # new code
    doublets_used = singlets_pairs_df.head(num_doublets)
    doublet_id = pd.DataFrame(doublets_used[0] + "--" + doublets_used[1])
    singlets_used_for_doublets = [item for pair in zip(doublets_used[0], doublets_used[1]) for item in pair]
    singlets_to_keep_in_seurat = [x for x in singlets_in_seurat if x not in singlets_used_for_doublets]
    if len(singlets_to_keep_in_seurat) > (num_doublets * (1 - doublets_pct) / doublets_pct):
        singlets_trimmed = random.sample(singlets_to_keep_in_seurat,
                                         int(num_doublets * (1 - doublets_pct) / doublets_pct))
    else:
        num_doublets = int(len(singlets_in_seurat) * doublets_pct / (
                    1 - doublets_pct))  # calculate the number of doublets, no assumptions about removal being made
        avg_assay = avg_assay_ori[:num_doublets, :]
        doublets_used = singlets_pairs_df.head(num_doublets)
        singlets_used_for_doublets = [item for pair in zip(doublets_used[0], doublets_used[1]) for item in pair]
        singlets_to_keep_in_seurat = [x for x in singlets_in_seurat if x not in singlets_used_for_doublets]
        num_doublets_adjusted = int((len(singlets_to_keep_in_seurat) * doublets_pct) / (
                    1 - doublets_pct))  # adjust the number of doublets based on the number of singlets removed
        avg_assay = avg_assay[:num_doublets_adjusted,
                    :]  # subset the avg_assay_new to be as long as the adjusted doublet count
        doublets_used = singlets_pairs_df.head(num_doublets_adjusted)
        doublet_id = pd.DataFrame(doublets_used[0] + "--" + doublets_used[1])
        singlets_trimmed = singlets_to_keep_in_seurat  # no more singlets removed after the ones used to make doublets are removed

    singlets_trimmed_idx = [singlets_in_seurat.index(x) for x in singlets_trimmed]
    # switch the order to be the order of singlets trimmed
    counts_matrix = counts_matrix[singlets_trimmed_idx, :]

    doublet_labels = ["doublet" for i in range(avg_assay.shape[0])]
    singlet_labels = ["singlet" for i in range(counts_matrix.shape[0])]
    doublet_id = doublet_id[0].tolist()
    doublet_labels.extend(singlet_labels)
    doublet_assay = vstack((avg_assay, counts_matrix))
    doublet_id.extend(singlets_trimmed)
    estimated_doublet_rate = avg_assay.shape[0] / (avg_assay.shape[0] + counts_matrix.shape[0])
    print(estimated_doublet_rate)
    #     print("Total doublets: {}".format(str(avg_assay.shape[0])))
    #     print("Total dataset: {}".format(str(doublet_assay.shape[0])))
    # random initiation
    cur_random_state = random.randint(1, int(1e5))
    print("INFO: Random stat is {}".format(str(cur_random_state)))
    cur_scrub = scr.Scrublet(doublet_assay, expected_doublet_rate=exp_dbl, random_state=cur_random_state)
    # Unlike Seurat, we must return the labels separately instead of as a part of metadata
    return cur_scrub, doublet_labels, doublet_id


def run_scrublet(data_dir, exp_dbl, act_dbl, save_dir):
    # ValueError: n_components=30 must be between 1 and min(n_samples, n_features)=9 with svd_solver='arpack'
    cur_scrub, doublet_labels, barcodes = load_fatemap_into_scrublet(data_dir,
                                                           exp_dbl=exp_dbl, doublets_pct=act_dbl)

    doublet_scores, predicted_doublets = cur_scrub.scrub_doublets()
    df = pd.DataFrame({'score': doublet_scores, 'label': predicted_doublets, 'barcode': barcodes})
    prefix1 = os.path.basename(data_dir)
    parts = data_dir.split('/')
    prefix2 = parts[-3]
    outfile = f"{save_dir}/{prefix2}___{prefix1}___exp_{exp_dbl}__act_{act_dbl}.csv"
    df.to_csv(outfile, index = False)
    fpr, tpr, threshold = metrics.roc_curve(y_true=doublet_labels, y_score=doublet_scores, pos_label="doublet")
    roc_auc = metrics.auc(fpr, tpr)
    precision, recall, thresholds = metrics.precision_recall_curve(doublet_labels, doublet_scores, pos_label="doublet")
    auc_precision_recall = metrics.auc(recall, precision)
    out_df = pd.DataFrame({
        "roc": [roc_auc],
        "pr": [auc_precision_recall]
    })
    return out_df


def run_bchmk(data_dir, data_sample_ID, num_to_test, exp_dbl, act_dbl, save_dir):
    df_ls = []
    for idx in range(num_to_test):
        cur_res = run_scrublet(data_dir, exp_dbl, act_dbl, save_dir)
        df_ls.append(cur_res)

    all_stats = pd.concat(df_ls)
    all_stats['ID'] = data_sample_ID
    return all_stats


# need to reload scrublet object everytime it is ran because it does not reset
def run_scrublet_on_dataset(data_dir, METHOD_ID, exp_dbl, act_dbl, num_to_test=5,
                            save_root="/projects/b1042/GoyalLab/zzhang/doublet_objects/scrublet"):
    print("Working on dataset: {}".format(data_dir))
    dir_10X = data_dir + "/10X"
    dataset_id = data_dir.split("/")[-1]
    sample_dirs = [dir_10X + "/" + x for x in os.listdir(dir_10X)]
    out_df = pd.DataFrame(columns=["roc", "pr"])
    for cur_sample_dir in sample_dirs:
        cur_sample_id = cur_sample_dir.split("/")[-1]
        cur_data_sample_id = "{}_{}".format(dataset_id, cur_sample_id)
        act_dbl_str = str(act_dbl)
        cur_save_dir = f"{save_root}/act__{act_dbl_str}"
        if not os.path.exists(cur_save_dir):
            os.makedirs(cur_save_dir)
        cur_df = run_bchmk(cur_sample_dir, cur_data_sample_id, num_to_test, exp_dbl, act_dbl,
                           cur_save_dir)
        out_df = pd.concat([out_df, cur_df])
    return out_df
