import numpy as np
merged_dict = {}


file_1 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_0_1000.npy",
        allow_pickle=True).item()
print(type(file_1))
file_2 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_1000_2000.npy",
        allow_pickle=True).item()

file_3 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_2k_2k5.npy",
        allow_pickle=True).item()

file_4 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_2k5_3k.npy",
        allow_pickle=True).item()

file_5 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_3k_3k5.npy",
        allow_pickle=True).item()
file_6 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_3k5_4k.npy",
        allow_pickle=True).item()
file_7 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_4k_4k5.npy",
        allow_pickle=True).item()

file_8 = np.load(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/attention_4k5_4804.npy",
        allow_pickle=True).item()

merged_dict.update(file_1)
merged_dict.update(file_2)
merged_dict.update(file_3)
merged_dict.update(file_4)
merged_dict.update(file_5)
merged_dict.update(file_6)
merged_dict.update(file_7)
merged_dict.update(file_8)


np.save("/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/attention_0_4804", merged_dict)


