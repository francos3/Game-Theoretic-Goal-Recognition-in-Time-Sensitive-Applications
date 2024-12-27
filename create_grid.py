filename="Grid-"
n=3
import random
def create_grid_graph(n):
    global filename
    filename+=str(n)+".csv"
    # Initialize an empty list of edges
    edges = []
    cost=1

    # Loop over each cell in the grid
    for i in range(n):
        for j in range(n):
            #cost=random.randint(1,100)
            # Calculate the node number for this cell
            node = i * n + j

            # Check the cell above if it exists
            if i > 0:
                edges.append((node, (i-1)*n + j,cost))
            # Check the cell below if it exists
            if i < n-1:
                edges.append((node, (i+1)*n + j,cost))
            # Check the cell to the left if it exists
            if j > 0:
                edges.append((node, i*n + j-1,cost))
            # Check the cell to the right if it exists
            if j < n-1:
                edges.append((node, i*n + j+1,cost))
    return edges

def save_graph_to_file(edges, filename):
    with open(filename, 'w') as file:
        for edge in edges:
            print(edge)
            file.write(','.join(map(str, edge)) + '\n')

# Example usage:
edges = create_grid_graph(n)
save_graph_to_file(edges, filename)
