import numpy as np
from Bio import SeqIO
import pandas as pd
import random

BINS = 100
THRESHOLD = 0.003

if __name__ == "__main__":

    dict_out1 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_0_4804.npy",
        allow_pickle=True).item()

    pos_records = list(SeqIO.parse(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/m6A_3UTR_ind_pos.fa",
        "fasta"))
    neg_records = list(SeqIO.parse(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/m6A_3UTR_ind_neg_20_window.fa",
        "fasta"))

    # counter = 0
    # pos_attention_scores = []
    # pos_label = []
    # neg_attention_scores = []
    # neg_label = []
    #
    # counter_10 = 0
    # counter_20 = 0
    # counter_30 = 0
    # counter_40 = 0
    # counter_50 = 0
    # counter_60 = 0
    # counter_70 = 0
    # counter_80 = 0
    # counter_90 = 0
    # counter_100 = 0
    # counter_10_origin = 0
    # counter_20_origin = 0
    # counter_30_origin = 0
    # counter_40_origin = 0
    # counter_50_origin = 0
    # counter_60_origin = 0
    # counter_70_origin = 0
    # counter_80_origin = 0
    # counter_90_origin = 0
    # counter_100_origin = 0

    # counter_all = 0
    # counter_threshold = 0
    all_sub = []
    all_sub_2 = []
    result_dict_all = {}
    result_dict_threshold = {}
    # result_all_actual_data = []
    # result_threshold_actual_data = []
    for a in range(BINS):
        result_dict_all["bin{bin_num}".format(bin_num=a + 1)] = 0
        result_dict_threshold["bin{bin_num}".format(bin_num=a + 1)] = 0

    for i in range(len(pos_records)):
        print(i)
        print(len(pos_records))
        print(pos_records[i].description.split(" "))
        pos_index_pre = pos_records[i].description.split(" ")[-1]
        pos_index = pos_index_pre.split("/")
        # neg_index_pre = neg_records[i].description.split(" ")[-1]
        # neg_index = neg_index_pre.split("/")
        # print("neg_index: ", neg_index)
        seq_len = len(str(pos_records[i].seq))
        all_index = pos_index  # + neg_index
        for ind in all_index:
            for bin_own in range(BINS):
                if int(ind) <= seq_len * ((bin_own + 1) / BINS):
                    result_dict_all[
                        "bin{bin_num}".format(bin_num=bin_own + 1)] = \
                        result_dict_all[
                            "bin{bin_num}".format(bin_num=bin_own + 1)] + 1
                    if THRESHOLD <= dict_out1[pos_records[i].id][int(ind)]:
                        result_dict_threshold[
                            "bin{bin_num}".format(bin_num=bin_own + 1)] = \
                            result_dict_threshold[
                                "bin{bin_num}".format(bin_num=bin_own + 1)] + 1

                        all_sub.append(random.uniform(((1 / BINS) * bin_own) + 0.01,
                                                      (1 / BINS) * (bin_own + 1)))
                    break

            # if int(ind) <= seq_len * 0.1:
            #     counter_10 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_10_origin += 1
            #         all_sub.append(random.uniform(0,0.1))
            # elif int(ind) <= seq_len * 0.2:
            #     counter_20 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_20_origin += 1
            #         all_sub.append(random.uniform(0.11,0.2))
            # elif int(ind) <= seq_len * 0.3:
            #     counter_30 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_30_origin += 1
            #         all_sub.append(random.uniform(0.21,0.3))
            # elif int(ind) <= seq_len * 0.4:
            #     counter_40 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_40_origin += 1
            #         all_sub.append(random.uniform(0.31,0.4))
            # elif int(ind) <= seq_len * 0.5:
            #     counter_50 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_50_origin += 1
            #         all_sub.append(random.uniform(0.41,0.5))
            # elif int(ind) <= seq_len * 0.6:
            #     counter_60 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_60_origin += 1
            #         all_sub.append(random.uniform(0.51,0.6))
            # elif int(ind) <= seq_len * 0.7:
            #     counter_70 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_70_origin += 1
            #         all_sub.append(random.uniform(0.61,0.7))
            # elif int(ind) <= seq_len * 0.8:
            #     counter_80 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_80_origin += 1
            #         all_sub.append(random.uniform(0.71,0.8))
            # elif int(ind) <= seq_len * 0.9:
            #     counter_90 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_90_origin += 1
            #         all_sub.append(random.uniform(0.81,0.9))
            # elif int(ind) <= seq_len:
            #     counter_100 += 1
            #     if 0.003 <= dict_out1[pos_records[i].id][int(ind)]:
            #         counter_100_origin += 1
            #         all_sub.append(random.uniform(0.91,0.99))

    # result = []
    # result.append(counter_10)
    # result.append(counter_20)
    # result.append(counter_30)
    # result.append(counter_40)
    # result.append(counter_50)
    # result.append(counter_60)
    # result.append(counter_70)
    # result.append(counter_80)
    # result.append(counter_90)
    # result.append(counter_100)
    index = []
    for m in range(BINS):
        index.append((1/BINS)*m)

    result = []
    for b in range(BINS):
        result.append(result_dict_all[
                        "bin{bin_num}".format(bin_num=b + 1)])

    dataframe_normal = pd.DataFrame({"Frequency": result, "Relative_distance": index})
    dataframe_normal.to_csv(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_all.csv",
        index=False)

    result2 = []
    for c in range(BINS):
        result2.append(result_dict_threshold[
                        "bin{bin_num}".format(bin_num=c + 1)])
    dataframe_pos = pd.DataFrame(
        {"Frequency": result2, "Relative_distance": index})
    dataframe_pos.to_csv(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_threshold.csv",
        index=False)
    final_result = []
    for i in range(BINS):
        final_result.append(result2[i] / result[i])
    print(final_result)
    dataframe_final = pd.DataFrame(
        {"Frequency": final_result, "Relative_distance": index})
    dataframe_final.to_csv(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_tim_pos_final_normalized.csv",
        index=False)

    # dataframe_all_sub = pd.DataFrame(
    #     {'Relative_Distance': all_sub})
    # dataframe_all_sub.to_csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/experiment_tim/experiment_generated_all_points_threshold.csv")
