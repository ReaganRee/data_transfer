import numpy as np
from Bio import SeqIO
from tqdm import tqdm

if __name__ == "__main__":

    dict_out1 = np.load(
        "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/attention_4k5_4804.npy",
        allow_pickle=True).item()

    print("length:  ",dict_out1)

    pos_records = list(SeqIO.parse(
        "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/m6A_3UTR_ind_pos.fa",
        "fasta"))
    neg_records = list(SeqIO.parse(
        "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/m6A_3UTR_ind_neg_full.fa",
        "fasta"))
    counter = 0
    pos_attention_scores = []
    pos_label = []
    neg_attention_scores = []
    neg_label = []
    for i in range(4500,4804):
        print(i)
        print(len(pos_records))
        print(pos_records[i].description.split(" "))
        pos_index_pre = pos_records[i].description.split(" ")[-1]
        pos_index = pos_index_pre.split("/")
        neg_index_pre = neg_records[i].description.split(" ")[-1]
        neg_index = neg_index_pre.split("/")
        print("neg_index: ", neg_index)
        seq_len = len(str(pos_records[i].seq))
        for ind in pos_index:
            #print(pos_records[i].id)
            #print(int(ind) + 1)
            #if int(pos_index[n]) + 1 < seq_len and int(neg_index[n]) + 1 < seq_len:
            pos_score = dict_out1[pos_records[i].id][int(ind)]
            pos_attention_scores.append(pos_score * 100)
            pos_label.append("pos")

        for neg_ind in neg_index:
            print(neg_ind)
            print("seq_len: ", seq_len)
            print("attention_len: ", len(dict_out1[neg_records[i].id]))
            neg_score = dict_out1[neg_records[i].id][int(neg_ind)]
            neg_attention_scores.append(neg_score *100)
            neg_label.append("neg")
            # if pos_score > neg_score:
            #     counter += 1
                # print(pos_score - neg_score)

    import pandas as pd
    all_label = pos_label + neg_label
    all_scores = pos_attention_scores + neg_attention_scores
    dataframe_pos = pd.DataFrame({'pos_or_neg': all_label, "attention_score": all_scores})
    dataframe_pos.to_csv("/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/4k5_4804_attention.csv", index=False)

    # dataframe_neg = pd.DataFrame({'neg_attention_scores': neg_attention_scores})
    # dataframe_neg.to_csv(
    #     "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/neg_0_50_attention.csv",
    #     index=False)

    # print("reagan: ",counter)
    # print("total_counter: ", total_counter)

    # print(dict_out1['ENST00000400830.4|ENSG00000197785.14|OTTHUMG00000000575.7|OTTHUMT00000100388.4|ATAD3A-204|ATAD3A|830|CDS:1-760|UTR3:761-830|'][40])
    # print(dict_out1[
    #           'ENST00000400830.4|ENSG00000197785.14|OTTHUMG00000000575.7|OTTHUMT00000100388.4|ATAD3A-204|ATAD3A|830|CDS:1-760|UTR3:761-830|'][16])

    # print(dict_out[
    #           'ENST00000308647.8|ENSG00000160072.20|OTTHUMG00000000577.6|-|ATAD3B-201|ATAD3B|2114|UTR5:1-101|CDS:102-995|UTR3:996-2114|'][
    #           626])
    # dict_out2 = np.load(
    #     "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/test_neg.npy",
    #     allow_pickle=True).item()
    # helper = []
    # for i in range(len(dict_out2[
    #           'ENST00000400830.4|ENSG00000197785.14|OTTHUMG00000000575.7|OTTHUMT00000100388.4|ATAD3A-204|ATAD3A|830|CDS:1-760|UTR3:761-830|'])):
    #     result = dict_out2[
    #         'ENST00000400830.4|ENSG00000197785.14|OTTHUMG00000000575.7|OTTHUMT00000100388.4|ATAD3A-204|ATAD3A|830|CDS:1-760|UTR3:761-830|'][i]- dict_out1['ENST00000400830.4|ENSG00000197785.14|OTTHUMG00000000575.7|OTTHUMT00000100388.4|ATAD3A-204|ATAD3A|830|CDS:1-760|UTR3:761-830|'][i]
    #     helper.append(result)
    # print(helper)
