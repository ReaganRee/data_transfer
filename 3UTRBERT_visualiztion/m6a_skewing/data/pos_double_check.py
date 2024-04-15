# from Bio import SeqIO
# import numpy as np
# # pos_records = list(SeqIO.parse(
# #         "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/m6A_3UTR_ind_pos.fa",
# #         "fasta"))
# # result = []
# # for seq_record in pos_records:
# #     pos_index_pre = seq_record.description.split(" ")[-1]
# #     pos_index = pos_index_pre.split("/")
# #     seq_len = len(str(seq_record.seq))
# #     for ind in pos_index:
# #         if int(ind) <= 20 or int(ind) >= seq_len - 20:
# #             result.append(seq_record.description)
# #
# # print(result)
#
# dict_out1 = np.load(
#         "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_0_4804.npy",
#         allow_pickle=True).item()
# #print(dict_out1)
# print(len(dict_out1["ENST00000652369.1|ENSG00000188157.15|OTTHUMG00000040778.7|OTTHUMT000"
#            "00097991.3|AGRN-210|AGRN|7408|UTR5:1-450|CDS:451-6273|UTR3:6274-7408|"]))
import numpy as np
A= []#np.array([[1,2,3,4], [5,6,7,8],[10,11,12,13]])
B= []#np.reshape(A, (3,4,1))
a= np.array([1,2,3,4])
A.append(a[:, np.newaxis])
A.append(a[:, np.newaxis])
print(np.array(A).shape)
print(A)
print(B.shape)
print(B)
