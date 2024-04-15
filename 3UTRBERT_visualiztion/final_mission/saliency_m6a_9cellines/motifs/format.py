with open("/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/motifs/MOLM13_motif_GACUUUC_1080.txt", "r") as file:
    row_list = file.readlines()
    with open("/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/sequences/MOLM13.txt", "w") as f2:
        for i in range(1,len(row_list), 2):
            f2.write(row_list[i])



