import csv
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import numpy as np


def process_csv(file_path, dictionary):
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        rows = [row for row in reader]
    #dictionary = {}
    for row in rows[1:]:
        string_checker = row[5]
        string_checker = string_checker.replace("A", "")
        string_checker = string_checker.replace("C", "")
        string_checker = string_checker.replace("G", "")
        string_checker = string_checker.replace("U", "")
        if string_checker == "" and row[4] == "Homo_sapiens":
            if row[2] + "_" + row[3] in dictionary:
                dictionary[row[2] + "_" + row[3]][0].append(row[5])
            else:
                dictionary[row[2] + "_" + row[3]] = [[row[5]], [], []]  #[[seq], [matrix_id], domain]
    return dictionary


def map_attract_db(dictionary):
    with open(
            "/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/ATtRACT/ATtRACT_db.txt",
            "r") as db:
        db_list = db.readlines()[1:]
        for key, value in dictionary.items():
            # print(key)
            for seq in value[0]:

                for motif in db_list:
                    row_list = motif.split("\t")
                    if row_list[1] == key.split("_")[-1] and seq in row_list[4]:
                        print("yes")
                        dictionary[key][1].append(row_list[-2]) # {gene_id: [seq, matrix_id]}
            # print(db_list)
    return dictionary


def map_domain(dictionary):
    with open(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/ATtRACT/domain.txt",
        "r") as domain:
        domain_list = domain.readlines()[1:-1]
        for key, value in dictionary.items():
            #print(key)
            #print(dictionary)
            for gene in domain_list:
                if gene.split("\t")[0] == key.split("_")[-1]:
                    dictionary[key][2].append(gene.split("\t")[1].strip())
        #print(dictionary)
    return dictionary


def map_matrix(dictionary, counter_dict):
    with open("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/ATtRACT/pwm.txt", "r") as file:
        pwm_list = file.readlines()
        #print(pwm_list)
    for key, value in dictionary.items():
        simplified_matrix_id = list(set(dictionary[key][1]))
        #counter_dict = {}
        for matrix_id in simplified_matrix_id:
            for pwm in pwm_list:
                if ">" in pwm and pwm.split("\t")[0] == ">" + matrix_id and 5< int(pwm.split("\t")[1]) < 8:
                    matrix = pwm_list[pwm_list.index(pwm) + 1:pwm_list.index(pwm) + 1 + int(pwm.split("\t")[1].strip())]
                    matrix_2d = []
                    for row in matrix:
                        row_list = row.split("\t")
                        matrix_2d.append(row_list)
                    domain = dictionary[key][2]
                    if "RRM" in domain:
                        filter_domain = "RRM"
                    elif "KH" in domain:
                        filter_domain = "KH"
                    elif "PUF" in domain:
                        filter_domain = "PUF"
                    elif "CSD" in domain:
                        filter_domain = "CSD"
                    else:
                        filter_domain = "Other"
                    new_matrix = np.zeros((4, int(pwm.split("\t")[1].strip())))
                    for i in range(4):
                        # row = []
                        for j in range(int(pwm.split("\t")[1].strip())):
                            new_matrix[i][j] = matrix_2d[j][i].strip()
                    #print(counter)
                    if key.split("_")[0] not in counter_dict:
                        counter_dict[key.split("_")[0]] = 1
                    else:
                        counter_dict[key.split("_")[0]] += 1
                    np.savetxt("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/motifs_PWM/{rbp}_{domain}_{num}.txt".format(
                                rbp=key.split("_")[0], domain=filter_domain, num=counter_dict[key.split("_")[0]]), new_matrix, fmt='%.2e')
                    #counter += 1
    return counter_dict




if __name__ == "__main__":
    counter_dict = {}
    dictionary = {}
    for file in os.listdir("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/comparision_from_ATtRACT"):
        print(file)

        dictionary = process_csv("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/comparision_from_ATtRACT/" + file, dictionary)
        #print(dictionary1)
        dictionary = map_attract_db(dictionary)
        dictionary = map_domain(dictionary)
        #print(dictionary3)
    counter_dict = map_matrix(dictionary, counter_dict)
