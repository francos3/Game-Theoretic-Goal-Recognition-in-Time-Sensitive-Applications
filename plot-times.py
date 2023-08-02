import re
import glob
import statistics
import csv
from collections import defaultdict

def extract_time(file_paths):
    time_categories = defaultdict(list)
    pattern = re.compile(r'log_prefix_N(?P<N>\d+)_R(?P<R>\d+)_S.*\.txt')
    time_pattern = re.compile(r'main,time:,(?P<time>\d+\.\d+)')

    for file_path in file_paths:
        try:
            match = pattern.search(file_path)
            if match:
                N = int(match.group('N'))
                R = int(match.group('R'))
                category = (N, R)

                with open(file_path, 'r') as file:
                    for line in file:
                        time_match = time_pattern.search(line)
                        if time_match:
                            time_value = float(time_match.group('time'))
                            time_categories[category].append(time_value)

        except FileNotFoundError:
            print(f"File {file_path} not found.")
            continue

    return time_categories

# Folder path
folder_path = '/home/ugac002/Dropbox/GO_2023/test_experiments/'

# Get file paths matching the given format
file_paths = glob.glob(f'{folder_path}log_prefix_N*_R*_S*.txt')

time_categories = extract_time(file_paths)

# Write to CSV file
with open('times_statistics.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['N', 'R', 'Average', 'Standard Deviation'])

    for category, times in time_categories.items():
        N, R = category
        if len(times)<2:
            continue
        average_time = round(statistics.mean(times),2)
        print("times:",times)
        standard_deviation_time = round(statistics.stdev(times))
        writer.writerow([N, R, average_time, standard_deviation_time])

import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('times_statistics.csv')

# Pivot the data to get R values as columns and N values as indices
pivot_data = data.pivot(index='N', columns='R', values='Average')

# Reset index to make N a column
pivot_data.reset_index(inplace=True)

# Melt the pivot data back to long format, keeping NaN for missing R values
melted_data = pivot_data.melt(id_vars=['N'], value_vars=pivot_data.columns[1:],
                              var_name='R', value_name='Average')

# Group by the N value
grouped_data = melted_data.groupby('N')

# Create a plot for each N value, adding a marker for each populated data point
for name, group in grouped_data:
    plt.plot(group['R'], group['Average'], label=f'N={name}', marker='o')

# Set x-ticks to the unique R values
plt.xticks(data['R'].unique())

# Add labels, title, and legend
plt.xlabel('Destinations')
plt.ylabel('Average time (secs)')
plt.title('Runtime as a function of grid size and |D|')
plt.legend()

# Add a grid
plt.grid(True)

# Show the plot
plt.show()