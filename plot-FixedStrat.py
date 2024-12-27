import re
import glob
import statistics
import csv
from collections import defaultdict
from pprint import pprint  # Import pprint for pretty printing dictionaries
from matplotlib.ticker import MaxNLocator


skip_extract=False

def extract_data(file_paths):
    if skip_extract:
        return
    data_categories = defaultdict(lambda: defaultdict(list))
    #pattern = re.compile(r'log_prefix_N(?P<N>\d+)_R(?P<R>\d+)_F(?P<F>\d+)_P(?P<P>\d+)_S.*\.txt')
    pattern = re.compile(r'log_prefix_NShanghai_R(?P<R>\d+)_F(?P<F>\d+)_P(?P<P>\d+)_S.*\.txt')
    #data_pattern = re.compile(r'avg_iterative_reward:,(?P<avg_iterative_reward>-?\d+(\.\d+)?([eE][-+]?\d+)?),stdev:,(?P<stdev>-?\d+(\.\d+)?([eE][-+]?\d+)?),CV:, ?(?P<CV>-?\d+(\.\d+)?([eE][-+]?\d+)?),entries:(?P<entries>\d+)')
    data_pattern = re.compile(r'avg_iterative_reward:,(?P<avg_iterative_reward>-?\d+(\.\d+)?([eE][-+]?\d+)?)')
    time_mem_pattern = re.compile(r'FINISHED CPU (?P<cpu_time>\d+\.\d+) MEM \d+ MAXMEM (?P<max_mem>\d+)')


    for file_path in file_paths:
        try:
            match = pattern.search(file_path)
            if match:
                #N = int(match.group('N'))
                R = int(match.group('R'))
                F = int(match.group('F'))
                P = int(match.group('P'))
                #category = (N, F, R, P)
                category = (F, R, P)
                
                with open(file_path, 'r') as file:
                    #print("working on file:",file)
                    avg_iterative_reward_found=False
                    for line in file:
                        data_match = data_pattern.search(line)
                        if data_match:
                            avg_iterative_reward_found=True
                            avg_iterative_reward = float(data_match.group('avg_iterative_reward'))

                            data_categories[category]['avg_iterative_reward'].append(avg_iterative_reward)
                        
                        time_mem_match = time_mem_pattern.search(line)
                        if time_mem_match:
                            cpu_time = float(time_mem_match.group('cpu_time'))
                            max_mem = int(time_mem_match.group('max_mem'))

                    if avg_iterative_reward_found:
                        data_categories[category]['cpu_time'].append(cpu_time)
                        data_categories[category]['max_mem'].append(max_mem/1000.0)
                    else:
                        print("avg_iterative_reward not found for file:",file)
        except FileNotFoundError:
            print(f"File {file_path} not found.")
            continue

    return data_categories
# Folder path
folder_path = '/home/ugac002/Dropbox/GO_2023/experiments_FictitiousPlay_Rationalized/'
# Get file paths matching the given format
#file_paths = glob.glob(f'{folder_path}log_prefix_N*_R*_S*.txt')
file_paths = glob.glob(f'{folder_path}log_prefix_N*_R*_S*.txt')
data_categories = extract_data(file_paths)

#pprint(dict(data_categories))

# Update the CSV writer section to include 'F' value
if not skip_extract:
    #sorted_categories = sorted(data_categories.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))  # Sort by N, F, R and P.
    sorted_categories = sorted(data_categories.keys(), key=lambda x: (x[0], x[1], x[2]))  # Sort by N, F, R and P.
    with open('data_statistics.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow(['N', 'F', 'R', 'P', 'Average Reward', 'Avg Time', 'Avg Mem'])
        writer.writerow(['F', 'R', 'P', 'Average Reward', 'Avg Time', 'Avg Mem'])

        for category in sorted_categories:
            values = data_categories[category]
            #N, F, R, P = category
            F, R, P = category


            values = data_categories[category]
            if all(len(v) < 2 for v in values.values()):
                #N,F, R, P = category
                F, R, P = category
                #print("\t skipping N:,",N,"F:,",F,"R:,",R, "P:",P)
                print("\t skipping F:,",F,"R:,",R, "P:",P)
                continue

            #print("N:,", N, "F:,",F,"R:,", R, "P:,",P, "data:,", len(values['avg_iterative_reward']))
            print("F:,",F,"R:,", R, "P:,",P, "data:,", len(values['avg_iterative_reward']))
            average_mean_savings = round(statistics.mean(values['avg_iterative_reward']), 2)
            average_time = round(statistics.mean(values['cpu_time']), 2)
            average_mem = round(statistics.mean(values['max_mem']), 2)
            #writer.writerow([N, F, R, P,average_mean_savings, average_time, average_mem])
            writer.writerow([F, R, P,average_mean_savings, average_time, average_mem])

#Plotting
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Fully use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{nccmath}'

# Read the CSV file
data = pd.read_csv('data_statistics.csv')

# Filter data based on unique N values
#unique_N_values = data['N'].unique()

# Create a new column 'F,R' by converting 'F' and 'R' columns to string and concatenating them
data['F,R,P'] = data['F'].astype(str) + ',' + data['R'].astype(str) + data['P'].astype(str) 

for (p_value, f_value), subset in data.groupby(['P', 'F']):
    # Plotting the subset
    plt.plot(subset['R'], subset['Average Reward'], marker='o', label=f'P={p_value}, K={f_value}')
    # Annotating each point
    for idx, row in subset.iterrows():
        x_value, y_value = row['R'], row['Average Reward']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(0, -10), ha='right', fontsize=7)

plt.xlabel('Destinations')
plt.ylabel('Average Reward')
plt.title('Fixed Stratedy, Average reward per |D|, Saturation Parameter K and Radius P.\n30x30 orthogonal graph with randomized origin and destinations.\n'+r'$q(d_{i})=min(1.0,\frac{1}{MaxCost(d_{i})/k})$')
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
#plt.legend(loc='upper left', bbox_to_anchor=(0.70, 0.98), borderaxespad=0.,ncol=2)  # Position the legend to the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(0.70, 0.98), borderaxespad=0.)  # Position the legend to the right of the plot
plt.grid(True)
plt.savefig('average_reward_per_FandP.png')
plt.show()

# Save the plot to a file, if needed



exit(1)

for N in unique_N_values:
    filtered_data = data[data['N'] == N]
   
    plt.figure(figsize=(12, 8))  # Create a new figure with a specific size
   #IndividualRewardPlot
    plt.plot(filtered_data['R'], filtered_data['Average Reward'], marker='o')
    plt.title('Average Reward, q=min(1.0,weight/F),1000 Paths per Destination')
    plt.xlabel('Destinations')
    plt.ylabel('Average Reward', fontweight='bold')  # Make y-label bold
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylim(bottom=min(filtered_data['Average Reward']) - 0.05, top=max(filtered_data['Average Reward']) + 0.05)  # Extend y-axis limits

    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Average Reward']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(0, -10), ha='right', fontsize=7)

    plt.savefig(f'AvgReward_for_N_{N}.png', dpi=300)  # Save the plot as a PNG with higher DPI
    plt.show()
    exit(1)
    
    plt.figure(figsize=(12, 8))  # Create a new figure with a specific size
    plt.suptitle('30x30 graph,F=2,q=min(1.0,weight/F) 1000 Paths per Destination', fontsize=16, fontweight='bold')
    


    # Plot for Average Reward
    plt.subplot(221)  # 2x2 grid, 1st subplot
    plt.plot(filtered_data['R'], filtered_data['Average Reward'], marker='o')
    plt.title(f'Average Reward,F=2,q=min(1.0,weight/F),1000 Paths per Destination')
    plt.xlabel('Destinations')
    plt.ylabel('Average Reward', fontweight='bold')  # Make y-label bold
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylim(bottom=min(filtered_data['Average Reward']) - 0.05, top=max(filtered_data['Average Reward']) + 0.05)  # Extend y-axis limits
    #plt.yticks(np.arange(min(filtered_data['Average Reward']), max(filtered_data['Average Reward']), 50.0))

    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Average Reward']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(0,-10), ha='right', fontsize=7)


    # Plot for CV
    plt.subplot(222)  # 2x2 grid, 2nd subplot
    plt.plot(filtered_data['R'], filtered_data['Average CV'], marker='o')
    plt.title(f'Average CV')
    plt.xlabel('Destinations')
    plt.ylabel('Average CV', fontweight='bold')  # Make y-label bold
    plt.yticks(np.arange(min(filtered_data['Average CV']), max(filtered_data['Average CV'])+0.01, 0.01))
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylim(bottom=min(filtered_data['Average CV']) - 0.001, top=max(filtered_data['Average CV']) + 0.006)  # Extend y-axis limits
    #plt.xlim(left=min(filtered_data['R']) - 1, right=max(filtered_data['R']) + 1)  # Extend x-axis limits
    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Average CV']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(5,5), ha='right',fontsize=7)


    # Plot for Distinct Cutsets
    plt.subplot(223)  # 2x2 grid, 3rd subplot
    plt.plot(filtered_data['R'], filtered_data['Distinct Cutsets'], marker='o')
    plt.title(f'Distinct Cutsets')
    plt.xlabel('Destinations')
    plt.ylabel('Distinct Cutsets', fontweight='bold')  # Make y-label bold
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylim(bottom=min(filtered_data['Distinct Cutsets']) - 2.0, top=max(filtered_data['Distinct Cutsets']) + 6.0)  # Extend y-axis limits
    #plt.ylim(bottom=min(filtered_data['Average Reward']) - 5, top=max(filtered_data['Average Reward']) + 5)  # Extend y-axis limits
    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Distinct Cutsets']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(5,5), ha='right',fontsize=7)
    
    # Plot for Last Cutset Size
    plt.subplot(224)  # 2x2 grid, 4th subplot
    plt.plot(filtered_data['R'], filtered_data['Last Cutset Size'], marker='o')
    plt.title(f'Last Cutset Size')
    plt.xlabel('Destinations')
    plt.ylabel('Last Cutset Size', fontweight='bold')  # Make y-label bold
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylim(bottom=min(filtered_data['Last Cutset Size']) - 0.5, top=max(filtered_data['Last Cutset Size']) + 0.5)  # Extend y-axis limits
    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Last Cutset Size']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(5,5), ha='right',fontsize=7)

    plt.tight_layout()
    #plt.ylim(bottom=min(filtered_data['Average Reward']) - 5, top=max(filtered_data['Average Reward']) + 5)  # Extend y-axis limits
    #plt.subplots_adjust(hspace=0.5, left=0.1, right=0.9, top=0.9, bottom=0.1)  # Adjust margins around subplots
    plt.subplots_adjust(hspace=0.5)  # Adjust margins around subplots

    plt.savefig(f'plot_for_N_{N}.png', dpi=300)  # Save the plot as a PNG with higher DPI
    plt.show()


#Now do time and memory
# New figure for Average Time and Average Memory
for N in unique_N_values:
    filtered_data = data[data['N'] == N]

    plt.figure(figsize=(12, 6))  # Create a new figure with a specific size
    plt.suptitle(f'Average Time and Memory for N = {N}, 1000 Paths per Destination.', fontsize=16, fontweight='bold')

    # Plot for Average Time
    plt.subplot(121)  # 1x2 grid, 1st subplot
    plt.plot(filtered_data['R'], filtered_data['Avg Time'], marker='o', color='b')
    plt.title(f'Average Time')
    plt.xlabel('Destinations')
    plt.ylabel('Average Time (s)', fontweight='bold')  # Make y-label bold
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Avg Time']
        plt.annotate(f'{y_value:.2f}', (x_value, y_value), textcoords="offset points", xytext=(0,10), ha='center', fontsize=7)

    # Plot for Average Memory
    plt.subplot(122)  # 1x2 grid, 2nd subplot
    plt.plot(filtered_data['R'], filtered_data['Avg Mem'], marker='o', color='r')
    plt.title(f'Average Memory')
    plt.xlabel('Destinations')
    plt.ylabel('Average Memory (MB)', fontweight='bold')  # Make y-label bold
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    for idx, row in filtered_data.iterrows():
        x_value, y_value = row['R'], row['Avg Mem']
        plt.annotate(f'{y_value}', (x_value, y_value), textcoords="offset points", xytext=(0,10), ha='center', fontsize=7)

    #plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # To make sure the suptitle fits
    plt.savefig(f'avg_time_mem_for_N_{N}.png', dpi=300)  # Save the plot as a PNG with higher DPI
    plt.show()