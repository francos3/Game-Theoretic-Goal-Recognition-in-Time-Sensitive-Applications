import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv('CombinedResultsN30.csv')

# Define line styles based on category keywords
line_styles = {
    'Fixed': 'solid',    # Solid line for "Fixed"
    'Reactive': 'dashed',  # Dashed line for "Reactive"
    'A Priori': 'dotted'  # Dotted line for "A Priori"
}

# Create a figure and plot
plt.figure(figsize=(12, 8))
for category, group in data.groupby('Category'):
    # Determine line style based on category name
    if 'Fixed' in category:
        linestyle = line_styles['Fixed']
    elif 'Reactive' in category:
        linestyle = line_styles['Reactive']
    elif 'A Priori' in category:
        linestyle = line_styles['A Priori']
    else:
        linestyle = '-'  # Default linestyle if no specific keyword found
    
    # Sort the group by 'Destinations' for coherent plotting
    group = group.sort_values(by='Destinations')
    plt.plot(group['Destinations'], group['Average Reward'], marker='o', linestyle=linestyle, label=category)

# Adding labels and title
plt.xlabel('Destinations')
plt.ylabel('Average Reward')
plt.title('Average Reward by Destinations for Each Category')
plt.legend(title='Category', loc='upper right', bbox_to_anchor=(1.25, 1))
plt.grid(True)

# Save the plot as a PNG file
plt.savefig('/path/to/save/average_reward_per_category.png', bbox_inches='tight')

# Show the plot
plt.show()

