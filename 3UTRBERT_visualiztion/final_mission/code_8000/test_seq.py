import os
from Bio import SeqIO
#print(len("abcde"))
for file in os.listdir("/Users/reagan/Desktop/3UTRBERT_visualiztion/final_mission/code_8000/5_fold_data_final"):
    if file.split(".")[-1] == "fasta":
        print(file)
        records_list = list(SeqIO.parse(
                "/Users/reagan/Desktop/3UTRBERT_visualiztion/final_mission/code_8000/5_fold_data_final/" + file,
                "fasta"))
        for seq_record in records_list:
            if len(str(seq_record.seq)) > 6000:
                print("YYYYYYYYYYYYY!")
