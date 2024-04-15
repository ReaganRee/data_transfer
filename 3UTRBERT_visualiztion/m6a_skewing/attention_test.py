import csv
import numpy as np
from scipy.stats import mannwhitneyu


pos_samples = []
neg_samples = []
count = 0
with open('/Users/reagan/Desktop/3UTRBERT_visualiztion/m6a_skewing/data/original/no_window/no_merged.csv') as f:
    for row in csv.reader(f, skipinitialspace=True):
        if row[0] == 'pos':
            count += 1
        if row[0] in ['pos', 'neg'] and row[1] != '' and float(row[1]) <= 0.5:
            if row[0] == 'pos':
                pos_samples.append(float(row[1]))
            elif row[0] == 'neg':
                neg_samples.append(float(row[1]))
        #print(row)


    positive_samples = np.array(pos_samples)
    negative_samples = np.array(neg_samples)
    pos_mean = np.mean(positive_samples)
    neg_mean = np.mean(negative_samples)
    stat, p_value = mannwhitneyu(positive_samples, negative_samples)

    print("U statistic:", stat)
    print("p-value:", p_value)
    print("pos mean:", pos_mean)
    print("neg mean:", neg_mean)
    print(count)
f.close()
