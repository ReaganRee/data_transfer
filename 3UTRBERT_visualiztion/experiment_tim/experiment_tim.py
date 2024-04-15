import numpy as np
from Bio import SeqIO
import pandas as pd
import random

if __name__ == "__main__":

    dict_out1 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_0_1000.npy",
        allow_pickle=True).item()

    pos_records = list(SeqIO.parse(
        "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/m6A_3UTR_ind_pos.fa",
        "fasta"))
    neg_records = list(SeqIO.parse(
        "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/m6A_3UTR_ind_neg_20_window.fa",
        "fasta"))

    counter = 0
    pos_attention_scores = []
    pos_label = []
    neg_attention_scores = []
    neg_label = []
    all_sub = []
    all_sub_normalized = []

    for i in range(len(pos_records)):
        print(i)
        print(len(pos_records))
        print(pos_records[i].description.split(" "))
        pos_index_pre = pos_records[i].description.split(" ")[-1]
        pos_index = pos_index_pre.split("/")
        neg_index_pre = neg_records[i].description.split(" ")[-1]
        neg_index = neg_index_pre.split("/")
        print("neg_index: ", neg_index)
        seq_len = len(str(pos_records[i].seq))
        all_index = pos_index # + neg_index
        for ind in all_index:
            if int(ind) <= seq_len * 0.1 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0,0.1))
            elif int(ind) <= seq_len * 0.2 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.11,0.2))
            elif int(ind) <= seq_len * 0.3 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.21,0.3))
            elif int(ind) <= seq_len * 0.4 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.31,0.4))
            elif int(ind) <= seq_len * 0.5 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.41,0.5))
            elif int(ind) <= seq_len * 0.6 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.51,0.6))
            elif int(ind) <= seq_len * 0.7 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.61,0.7))
            elif int(ind) <= seq_len * 0.8 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.71,0.8))
            elif int(ind) <= seq_len * 0.9 and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.81,0.9))
            elif int(ind) <= seq_len and 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
                all_sub.append(random.uniform(0.91,0.99))

        dataframe_pos = pd.DataFrame(
            {'points': all_sub})
        dataframe_pos.to_csv(
            "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/experiment_tim_pos.csv",
            index=False)



