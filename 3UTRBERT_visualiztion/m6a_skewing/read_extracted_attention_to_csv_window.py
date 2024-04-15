from statistics import mean

import numpy as np
from Bio import SeqIO
from tqdm import tqdm

if __name__ == "__main__":

    dict_out1 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_0_4804.npy",
        allow_pickle=True).item()

    #print("length:  ",dict_out1)

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
    for i in range(4500,4804):
        print(i)
        #print(len(pos_records))
        #print(pos_records[i].description.split(" "))
        pos_index_pre = pos_records[i].description.split(" ")[-1]
        pos_index = pos_index_pre.split("/")
        neg_index_pre = neg_records[i].description.split(" ")[-1]
        neg_index = neg_index_pre.split("/")
        #print("neg_index: ", neg_index)
        seq_len = len(str(pos_records[i].seq))
        for ind in pos_index:
            #print(pos_records[i].id)
            #print(int(ind) + 1)
            #if int(pos_index[n]) + 1 < seq_len and int(neg_index[n]) + 1 < seq_len:
            if int(ind) >= 20 and int(ind) < len(dict_out1[pos_records[i].id]) - 20:
                pos_score = mean(dict_out1[pos_records[i].id][int(ind)-20:int(ind) + 21])
                pos_attention_scores.append(pos_score * 100)
                pos_label.append("pos")

            else:
                pos_score = dict_out1[pos_records[i].id][int(ind)]
                pos_attention_scores.append(pos_score * 100)
                pos_label.append("pos")

        for neg_ind in neg_index:
            #print(neg_ind)
            # print("seq_len: ", seq_len)
            # print("attention_len: ", len(dict_out1[neg_records[i].id]))
            if int(neg_ind) >= 20 and int(neg_ind) < len(dict_out1[neg_records[i].id]) - 20:
                #print(dict_out1[neg_records[i].id][int(neg_ind)-5:int(neg_ind) + 6])
                neg_score = mean(dict_out1[neg_records[i].id][int(neg_ind)-20:int(neg_ind) + 21])
                neg_attention_scores.append(neg_score * 100)
                neg_label.append("neg")
            else:
                neg_score = dict_out1[neg_records[i].id][int(neg_ind)]
                neg_attention_scores.append(neg_score * 100)
                neg_label.append("neg")
            # if pos_score > neg_score:
            #     counter += 1
                # print(pos_score - neg_score)

    import pandas as pd
    all_label = pos_label + neg_label
    all_scores = pos_attention_scores + neg_attention_scores
    dataframe_pos = pd.DataFrame({'pos_or_neg': all_label, "attention_score": all_scores})
    dataframe_pos.to_csv("/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/4k5_4804_attention_window20_no_overlap.csv", index=False)

