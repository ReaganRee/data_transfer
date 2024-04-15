import numpy as np
from Bio import SeqIO
import pandas as pd
import random
BINS = 100
THRESHOLD = 0.003

dict_out1 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_0_4804.npy",
        allow_pickle=True).item()
print(dict_out1)
pos_records = list(SeqIO.parse(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/m6A_3UTR_ind_pos.fa",
    "fasta"))

attn_dict_all = {}
attn_dict_threshold = {}
for a in range(BINS):
    attn_dict_all["bin{bin_num}".format(bin_num=a + 1)] = 0
    attn_dict_threshold["bin{bin_num}".format(bin_num=a + 1)] = 0

for i in range(len(pos_records)):
    print(i)
    print(len(pos_records))
    print(pos_records[i].description.split(" "))
    pos_index_pre = pos_records[i].description.split(" ")[-1]
    pos_index = pos_index_pre.split("/")
    seq_len = len(str(pos_records[i].seq))
    all_index = pos_index
    for ind in all_index:
        for bin_own in range(BINS):
            if int(ind) <= seq_len * ((bin_own + 1) / BINS):
                attn_dict_all[
                    "bin{bin_num}".format(bin_num=bin_own + 1)] = \
                    attn_dict_all[
                        "bin{bin_num}".format(bin_num=bin_own + 1)] + dict_out1[pos_records[i].id][int(ind)]
                if THRESHOLD <= dict_out1[pos_records[i].id][int(ind)]:
                    attn_dict_threshold[
                        "bin{bin_num}".format(bin_num=bin_own + 1)] = \
                        attn_dict_threshold[
                            "bin{bin_num}".format(bin_num=bin_own + 1)] + dict_out1[pos_records[i].id][int(ind)]
                break
index = []
for m in range(BINS):
    index.append((1/BINS)*m)

result = []
for b in range(BINS):
    result.append(attn_dict_all[
                    "bin{bin_num}".format(bin_num=b + 1)])

dataframe_normal = pd.DataFrame({"attention_sum": result, "Relative_distance": index})
dataframe_normal.to_csv(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_all_attention.csv",
    index=False)

result2 = []
for c in range(BINS):
    result2.append(attn_dict_threshold[
                    "bin{bin_num}".format(bin_num=c + 1)])
dataframe_pos = pd.DataFrame(
    {"attention_sum": result2, "Relative_distance": index})
dataframe_pos.to_csv(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_threshold_attention.csv",
    index=False)
