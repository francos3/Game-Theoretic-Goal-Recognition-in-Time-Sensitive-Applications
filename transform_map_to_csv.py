import sys

def is_adjacent(grid, x, y, dx, dy):
    """Check if there is an adjacent node in the direction (dx, dy)."""
    new_x, new_y = x + dx, y + dy
    if 0 <= new_x < len(grid) and 0 <= new_y < len(grid[0]):
        return grid[new_x][new_y] == '.'
    return False

def parse_map_file(filename):
    """Parse the .map file starting from the fifth line and return a grid of nodes."""
    with open(filename, 'r') as file:
        lines = file.readlines()[4:]  # Skip the first four lines
        return [list(line.strip()) for line in lines]

def write_csv(filename, edges):
    """Write the edges to a CSV file."""
    with open(filename, 'w') as file:
        for edge in edges:
            file.write(f"{edge[0]},{edge[1]},1\n")

def main(map_filename):
    """Main function to read map file, process edges, and write to CSV."""
    grid = parse_map_file(map_filename)
    edges = set()
    directions = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
    
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] == '.':
                origin_node = f"{i},{j}"
                for dx, dy in directions:
                    if is_adjacent(grid, i, j, dx, dy):
                        destination_node = f"{i + dx},{j + dy}"
                        edge = (origin_node, destination_node)
                        edges.add(edge)  # Add each edge only once

    output_filename = map_filename.replace('.map', '.csv')
    write_csv(output_filename, edges)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py filename.map")
        sys.exit(1)

    map_filename = sys.argv[1]
    main(map_filename)

