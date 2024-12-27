import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import textwrap
import numpy as np



plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


# Load the CSV file
#data = pd.read_csv('CombinedResultsN30.csv')
data = pd.read_csv('CombinedResultsShanghai2.csv')

# Define line styles based on category keywords
line_styles = {
    'PB': 'dashed',    # Solid line for "Fixed"
    'NE': 'solid',  # Dashed line for "Best"
    'BL': 'dotted'  # Dotted line for "A Priori"
}



# Define specific colors for each 'K' value and 'A Priori'
category_colors = {
    't=1': 'blue',
    't=2': 'red',
    't=3': 'green',
    't=4': 'purple',
    'BL': 'gray'
}
# Function to wrap text for the legend
def wrap_legend_text(text, max_lines=4, width=20):
    wrapped_lines = textwrap.wrap(text, width=width)
    if len(wrapped_lines) > max_lines:
        wrapped_lines = wrapped_lines[:max_lines]
        wrapped_lines[-1] += '...'  # Indicate text was truncated
    return '\n'.join(wrapped_lines)

# Sort categories to ensure 'A Priori' is last
sorted_categories = sorted(data['Category'].unique(), key=lambda x: ('A Priori' in x, x))
# Create a figure and plot
plt.figure(figsize=(15, 8))

for category in sorted_categories:
    group = data[data['Category'] == category]
    # Determine line style based on category
    linestyle = line_styles.get('A Priori', '-')  # Default to solid if none of the keys match
    for key in line_styles:
        if key in category:
            linestyle = line_styles[key]
    
    # Determine color based on 'K' value or 'A Priori'
    color = 'gray'  # Default color if no match found
    for key in category_colors:
        if key in category:
            color = category_colors[key]
            break

    # Sort the group by 'Destinations' for coherent plotting
    group = group.sort_values('Destinations')

    # Wrap category text for legend
    #wrapped_category = wrap_legend_text(category, max_lines=4, width=20)

    plt.plot(group['Destinations'], group['Average Reward'], linestyle=linestyle, color=color, marker='o', label=category)

#Add dummy label to fix ordering in legends
plt.plot(2.5, np.zeros([1,3]), color='w', alpha=0, label=' ')

plt.ylim(0.1, 0.9)  # Increase upper limit by 0.1


# Adding labels and title
plt.xlabel(r'$|\mathcal{T}|$', fontsize=18, fontweight='bold')
plt.ylabel(r'Average $\Phi$',fontsize=18,fontweight='bold')
plt.tick_params(axis='both', labelsize=14)  # Increase font size to 14
#plt.title('Average '+r'$\phi(\mu,\rho)$' ' by Number of Destinations '+r'$|\mathcal{D}|$'+ ' for Each Category for a 30x30 orthogonal grid.\n'
#          'The Fixed Observer Strategy uses nodes up to half the optimal cost for each destination.\n'
#          ,fontsize=14,fontweight='bold')
#          #+
#          #'The Fixed Observer Strategy uses nodes up to half the optimal cost for each destination.\n'
#          #r'$q(w(\gamma^{k+})=min(1.0,\frac{w(\gamma^{k+}) \times t}{MaxWeight})$'+';'+r'$MaxWeight=\max\limits_{\gamma \in \mathcal{G}} (w(\gamma))$',fontsize=14,fontweight='bold'
#          #)
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

# Extract handles and labels after plotting
handles, labels = plt.gca().get_legend_handles_labels()

# Find index of 'BL'
index_of_BL = labels.index('BL')
bl_handle = handles.pop(index_of_BL)
bl_label = labels.pop(index_of_BL)

# Number of columns in the legend
num_columns = 3

# Calculate the position to insert 'BL' to start the third column
num_entries = len(handles) + 1  # Total entries including 'BL'
entries_per_column = (num_entries + num_columns - 1) // num_columns  # Ceil division
start_of_third_col = 2 * entries_per_column  # Index where third column starts

# Insert 'BL' at the start of the third column
handles.insert(start_of_third_col, bl_handle)
labels.insert(start_of_third_col, bl_label)

# Create a new legend with the reordered list
plt.legend(handles, labels, title='Categories', loc='upper left', bbox_to_anchor=(0.65, 1.02), ncol=3, fontsize=14)


plt.grid(True)

# Save the plot as a PNG file
plt.savefig('average_reward_per_category.png', bbox_inches='tight')

# Show the plot
plt.show()

