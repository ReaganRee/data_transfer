import csv

#def process_result(file_path):
























if __name__ == "__main__":
    import csv

    csv_reader = csv.reader(open("/Users/reagan/Desktop/final_mission/all_results/results/result_rnafm_22.csv"))
    for line in csv_reader:
        if line[1] == 'mean':
            mean_line = line[2:]
        print(line)
