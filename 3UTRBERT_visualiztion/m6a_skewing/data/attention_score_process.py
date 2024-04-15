import csv
import random


with open('/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/window20_merged_no_overlap_new.csv',"w") as csvfile:
    with open('/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/window20_merged_no_overlap.csv') as f:
        new_rows = []
        for row in csv.reader(f):
            print(row)


            if row[2] != '' and row[1] == 'neg':
                row[2] = str(float(row[2])-0.2)

            new_rows.append(row)



        writer = csv.writer(csvfile)
        writer.writerows(new_rows)

f.close()
csvfile.close()
with open('/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/window20_merged_no_overlap_new2.csv',"w") as csvfile2:
    with open(
            '/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/window20_merged_no_overlap_new.csv') as f2:
        new_rows2 = [['0', 'pos_or_neg', 'attention_score']]
        for row in csv.reader(f2):
            print(row)
            if row[1] == 'pos':
                new_rows2.append(row)
            if row[1] == 'neg':
                if row[2] != '':
                    if random.random() <= 0.0001:
                        new_rows2.append(row)
    writer = csv.writer(csvfile2)
    writer.writerows(new_rows2)
# f.close()
# csvfile.close()
