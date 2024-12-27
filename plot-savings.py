import re
import glob
import statistics
import csv
from collections import defaultdict

def extract_data(file_paths):
    data_categories = defaultdict(lambda: defaultdict(list))
    pattern = re.compile(r'log_prefix_N(?P<N>\d+)_R(?P<R>\d+)_S.*\.txt')
    #data_pattern = re.compile(r'mean_savings:,(?P<mean_savings>\d+\.\d+),std_dev:,(?P<std_dev>\d+\.\d+),initial_savings:,(?P<initial_savings>\d+\.\d+)')
    #data_pattern = re.compile(r'mean_savings:,(?P<mean_savings>\d+\.\d+),std_dev:,(?P<std_dev>\d+\.\d+),initial_savings:(?P<initial_savings>\d+\.\d+),dest_adj_savings:,(?P<dest_adj_savings>\d+\.\d+)')
    data_pattern = re.compile(r'mean_savings:,(?P<mean_savings>\d+\.\d+),std_dev:,(?P<std_dev>\d+\.\d+),initial_savings:(?P<initial_savings>\d+\.\d+),dest_adj_savings:,(?P<dest_adj_savings>\d+\.\d+),avg_target_choice:,(?P<avg_target_choice>\d+\.\d+)')


    for file_path in file_paths:
        try:
            match = pattern.search(file_path)
            if match:
                N = int(match.group('N'))
                R = int(match.group('R'))
                category = (N, R)
                #print("\t\tmatch found for N:",N,"R:",R)

                with open(file_path, 'r') as file:
                    for line in file:
                        data_match = data_pattern.search(line)
                        if data_match:
                            mean_savings = float(data_match.group('mean_savings'))
                            std_dev = float(data_match.group('std_dev'))
                            initial_savings = float(data_match.group('initial_savings'))
                            dest_adj_savings = float(data_match.group('dest_adj_savings'))
                            avg_target_choice = float(data_match.group('avg_target_choice'))
                            data_categories[category]['mean_savings'].append(mean_savings)
                            data_categories[category]['std_dev'].append(std_dev)
                            data_categories[category]['initial_savings'].append(initial_savings)
                            data_categories[category]['dest_adj_savings'].append(dest_adj_savings)
                            data_categories[category]['avg_target_choice'].append(avg_target_choice)
                            print("file:",file,"mean_savings:,",mean_savings,"avg_target_choice:,",avg_target_choice)

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
    writer.writerow(['N', 'R', 'Average Mean Savings', 'Average Std Dev', 'Average Initial Savings','Average DestAdj Savings','Average Target Choice'])

    for category, values in data_categories.items():
        #print("N:",N,"R:",R,"data:",len(values['mean_savings']))
        if all(len(v) < 2 for v in values.values()):
            print("\t skipping N:",N,"R:",R)
            continue

        N, R = category
        average_mean_savings = round(statistics.mean(values['mean_savings']), 2)
        average_std_dev = round(statistics.mean(values['std_dev']), 2)
        average_initial_savings = round(statistics.mean(values['initial_savings']), 2)
        average_dest_adj_savings = round(statistics.mean(values['dest_adj_savings']), 2)
        average_target_choice_savings = round(statistics.mean(values['avg_target_choice']), 2)
        writer.writerow([N, R, average_mean_savings, average_std_dev, average_initial_savings,average_dest_adj_savings,average_target_choice_savings])


############PLOTTING###################
import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('data_statistics.csv')

# Pivot the data to get R values as columns and N values as indices for both Average Mean Savings and Average Initial Savings
pivot_mean_savings = data.pivot(index='N', columns='R', values='Average Mean Savings')
pivot_initial_savings = data.pivot(index='N', columns='R', values='Average Target Choice')

# Reset index to make N a column
pivot_mean_savings.reset_index(inplace=True)
pivot_initial_savings.reset_index(inplace=True)

# Melt the pivot data back to long format, keeping NaN for missing R values
melted_mean_savings = pivot_mean_savings.melt(id_vars=['N'], value_vars=pivot_mean_savings.columns[1:],
                                              var_name='R', value_name='Average Mean Savings')
melted_initial_savings = pivot_initial_savings.melt(id_vars=['N'], value_vars=pivot_initial_savings.columns[1:],
                                                    var_name='R', value_name='Average Target Choice')

# Group by the N value for both melted dataframes
grouped_mean_savings = melted_mean_savings.groupby('N')
grouped_initial_savings = melted_initial_savings.groupby('N')

# Define some colors (you can modify this list to use your preferred colors)
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

# Create the figure with specified size
plt.figure(figsize=(8, 6))

# Create a plot for each N value for both Average Mean Savings and Average Target Choice
for (name1, group1), (name2, group2), color in zip(grouped_mean_savings, grouped_initial_savings, colors):
    plt.plot(group1['R'], group1['Average Mean Savings'], label=f'N={name1} Φ(C)', marker='o', color=color)
    plt.plot(group2['R'], group2['Average Target Choice'], label=f'N={name2} Φ(o)', marker='s', linestyle='dotted',color=color)

# Set y-axis to logarithmic scale
#plt.yscale('log')

# Set x-ticks to the unique R values
#print("data[R]:",data['R'])
plt.xticks(data['R'].unique())

# Add labels
plt.xlabel('Destinations')
plt.ylabel('Φ')

# Add legend above the plot in two rows
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.16), fontsize=9,ncol=4)

# Add a grid
plt.grid(True)

# Add title at the bottom
#plt.suptitle('Phi as a function of grid size and |D|', y=-0.1)
plt.title('Φ as a Function of Grid Size and Number of Destinations.', y=+1.00,fontweight="bold",fontsize=10)

# Adjust the layout
#plt.tight_layout()

# Show the plot
#plt.show()
plt.savefig("ComparativeRewards.png")
