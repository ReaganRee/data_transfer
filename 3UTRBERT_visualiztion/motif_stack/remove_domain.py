import os
import sys
folder = '/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/motifs_PWM'
import csv

headers = ["name", "domain"]
counter_dict = {}

with open("/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/domain_withGC.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(headers)
    for filename in os.listdir(folder):
        print(filename)
        print(counter_dict)
        if filename.split("_")[0] not in counter_dict:
            counter_dict[filename.split("_")[0]] = 1
            infilename = os.path.join(folder, filename)
            if not os.path.isfile(infilename):
                continue
            newname = infilename.replace(filename,
                                         filename.split("_")[0] + "_1" + ".pcm")
            writer.writerow([filename.split("_")[0] + "_1", filename.split("_")[1]])
            output = os.rename(infilename, newname)
        else:
            counter_dict[filename.split("_")[0]] += 1
            infilename = os.path.join(folder, filename)
            if not os.path.isfile(infilename):
                continue
            newname = infilename.replace(filename,
                                         filename.split("_")[0] + "_" + str(counter_dict[filename.split("_")[0]]) + ".pcm")
            writer.writerow([filename.split("_")[0] + "_" + str(counter_dict[filename.split("_")[0]]), filename.split("_")[1]])
            output = os.rename(infilename, newname)




