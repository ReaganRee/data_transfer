a = ['AKAP1_HepG2\n', 'BCLAF1_HepG2\n', 'DDX3X_HepG2\n', 'DDX3X_K562\n', 'DDX24_K562\n', 'FAM120A_K562\n', 'G3BP1_HepG2\n', 'GRWD1_HepG2\n', 'IGF2BP1_K562\n', 'LARP4_HepG2\n', 'LIN28B_K562\n', 'PABPC4_K562\n', 'PPIG_HepG2\n', 'PUM2_K562\n', 'RBM15_K562\n', 'RPS3_HepG2\n', 'SND1_HepG2\n', 'UCHL5_K562\n', 'UPF1_HepG2\n', 'UPF1_K562\n', 'YBX3_K562\n', 'ZNF622_K562']
b = []
for item in a:
    #item.strip()
    b.append(item.strip())
dictionary = dict.fromkeys(b, 0)
print(dictionary)
with open("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/for_ATtRACT.txt", "r") as f:
    pwm_list = f.readlines()
    for i in range(0, len(pwm_list), 6):
        rbp = str(pwm_list[i].split(" ")[1].split("_")[0]) + "_" + str(pwm_list[i].split(" ")[1].split("_")[1])
        dictionary[rbp] += 1

        with open("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/motifs_PWM/{name}_{num}_{motif}.pcm".format(name=rbp, num=str(dictionary[rbp]), motif=str(pwm_list[i].split(" ")[1].split("_")[2])), "w") as f2:
            f2.write(pwm_list[i+1])
            f2.write(pwm_list[i + 2])
            f2.write(pwm_list[i + 3])
            f2.write(pwm_list[i + 4])

