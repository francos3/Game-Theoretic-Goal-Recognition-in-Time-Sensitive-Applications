import re
import glob
import statistics
import csv
from collections import defaultdict

def extract_data(file_paths):
    data_categories = defaultdict(lambda: defaultdict(list))
    pattern = re.compile(r'log_prefix_N(?P<N>\d+)_R(?P<R>\d+)_S.*\.txt')
    data_pattern = re.compile(r'mean_savings:,(?P<mean_savings>\d+\.\d+),std_dev:,(?P<std_dev>\d+\.\d+),initial_savings:,(?P<initial_savings>\d+\.\d+)')

    for file_path in file_paths:
        try:
            match = pattern.search(file_path)
            if match:
                N = int(match.group('N'))
                R = int(match.group('R'))
                category = (N, R)

                with open(file_path, 'r') as file:
                    for line in file:
                        data_match = data_pattern.search(line)
                        if data_match:
                            mean_savings = float(data_match.group('mean_savings'))
                            std_dev = float(data_match.group('std_dev'))
                            initial_savings = float(data_match.group('initial_savings'))
                            data_categories[category]['mean_savings'].append(mean_savings)
                            data_categories[category]['std_dev'].append(std_dev)
                            data_categories[category]['initial_savings'].append(initial_savings)

        except FileNotFoundError:
            print(f"File {file_path} not found.")
            continue

    return data_categories

# Folder path
folder_path = '/home/ugac002/Dropbox/GO_2023/test_experiments/'

# Get file paths matching the given format
file_paths = glob.glob(f'{folder_path}log_prefix_N*_R*_S*.txt')

data_categories = extract_data(file_paths)

# Write to CSV file
with open('data_statistics.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['N', 'R', 'Average Mean Savings', 'Average Std Dev', 'Average Initial Savings'])

    for category, values in data_categories.items():
        if all(len(v) < 2 for v in values.values()):
            continue

        N, R = category
        average_mean_savings = round(statistics.mean(values['mean_savings']), 2)
        average_std_dev = round(statistics.mean(values['std_dev']), 2)
        average_initial_savings = round(statistics.mean(values['initial_savings']), 2)
        writer.writerow([N, R, average_mean_savings, average_std_dev, average_initial_savings])