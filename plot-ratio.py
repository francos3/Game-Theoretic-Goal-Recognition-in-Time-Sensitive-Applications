from collections import defaultdict
import glob
import re

# Folder path
folder_path = "/home/ugac002/Dropbox/GO_2023/test_experiments"

# Define a pattern to match the filenames
file_pattern = "log_prefix_N*_R*_S*.txt"

# Get all matching files
file_paths = glob.glob(f"{folder_path}/{file_pattern}")

# Dictionary to store percentages grouped by (N, R)
percentages_by_category = defaultdict(lambda: defaultdict(list))

# Iterate through the files and extract the data
for file_path in file_paths:
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Check if the file contains the "main,time:," line, and skip if not present
        if not any("main,time:" in line for line in lines):
            continue

        # Extract N and R values from the filename
        N_value = re.search(r'_N(\d+)_', file_path).group(1)
        R_value = re.search(r'_R(\d+)_', file_path).group(1)
        category = (N_value, R_value)

        # Extracting partial times for "Yen added"
        #print("filepath:",file_path)
        #yen_added_times = sum(float(re.search(r'paths in,(\d+.\d+)', line).group(1)) for line in lines if "Yen added" in line)
        yen_added_times = sum(float(re.search(r'paths in,(\d+\.?\d*)', line).group(1)) for line in lines if "Yen added" in line)


        # Extracting other specific time components
        link_prefix_time = float(re.search(r'linkPrefixToPaths,time:,(\d+.\d+)', "\n".join(lines)).group(1))
        #min_path_time = float(re.search(r'min_path_time,time:,(\d+.\d+)', "\n".join(lines)).group(1))
        calculate_q_recursive_time = float(re.search(r'calculate_q_recursive,time:,(\d+.\d+)', "\n".join(lines)).group(1))

        # Extracting total main time
        main_time = float(re.search(r'main,time:,(\d+.\d+)', "\n".join(lines)).group(1))

        # Calculate and store percentages separately
        percentages_by_category[category]['yen_added'].append(yen_added_times / main_time * 100)
        percentages_by_category[category]['link_prefix'].append(link_prefix_time / main_time * 100)
        #percentages_by_category[category]['min_path'].append(min_path_time / main_time * 100)
        percentages_by_category[category]['calculate_q_recursive'].append(calculate_q_recursive_time / main_time * 100)

# Calculate average percentages for each category
average_percentages_by_category = {category: {key: sum(percentages) / len(percentages) for key, percentages in time_parts.items()} for category, time_parts in percentages_by_category.items()}

#for category, average_percentages in average_percentages_by_category.items():
#    N_value, R_value = category
#    print(f"For N={N_value}, R={R_value}:")
#    for key, avg in average_percentages.items():
#        print(f"    {key}: {avg:.2f}%")


import csv

# Define the CSV filename
csv_filename = 'time_percentages.csv'

# Define the CSV headers
headers = ['N', 'R', 'Yen(%)', 'LinkPrefixToPaths(%)', 'CalculateQRecursive(%)']

# Open the CSV file for writing
with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write the headers to the CSV file
    writer.writerow(headers)

    # Iterate over the categories and average percentages, and write to the CSV file
    for category, average_percentages in average_percentages_by_category.items():
        N_value, R_value = category
        row = [
            N_value,
            R_value,
            round(average_percentages['yen_added'],2),
            round(average_percentages['link_prefix'],2),
            #round(average_percentages['min_path'],2),
            round(average_percentages['calculate_q_recursive'],2)
        ]
        writer.writerow(row)

print(f"Results saved to {csv_filename}")
import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('time_percentages.csv')

# List of unique N values
N_values = data['N'].unique()

# Define colors for different percentage categories
colors = {
    'Yen(%)': 'blue',
    'LinkPrefixToPaths(%)': 'green',
    'CalculateQRecursive(%)': 'red'
}

# Define line styles for different N categories
linestyles = ['-', '--', ':', '-.']
markers = ['o', 's', 'D', '^', 'v']

# Iterate through N values to plot separate lines for each N category
for idx, N in enumerate(N_values):
    subset = data[data['N'] == N]
    subset = subset.sort_values('R')

    # Define the linestyle based on the N category index
    linestyle = linestyles[idx % len(linestyles)]
    marker = markers[idx % len(markers)]


    # Plot each percentage category with the corresponding color and linestyle

    plt.plot(subset['R'], subset['Yen(%)'], color=colors['Yen(%)'], linestyle=linestyle, marker=marker, label=f'N={N} Yen(%)')
    plt.plot(subset['R'], subset['LinkPrefixToPaths(%)'], color=colors['LinkPrefixToPaths(%)'], linestyle=linestyle, marker=marker, label=f'N={N} LinkPrefixToPaths(%)')
    plt.plot(subset['R'], subset['CalculateQRecursive(%)'], color=colors['CalculateQRecursive(%)'], linestyle=linestyle, marker=marker, label=f'N={N} CalculateQRecursive(%)')

# Add a grid
plt.grid()

# Set labels and title
plt.xlabel("Destinations")
plt.ylabel("Percentage (%)")
#plt.title("Runtime as a function of grid size and |D|")

# Arrange the legend
plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.102), loc='lower left', ncol=2, mode="expand", borderaxespad=0,fontsize=8)

# Set x-ticks to only include populated R values
plt.xticks(data['R'].unique())

# Show the plot
plt.tight_layout()
plt.show()