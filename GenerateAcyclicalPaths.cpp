#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace std;

class Graph {
public:
    unordered_map<int, vector<pair<int, int>>> adj;

    void addEdge(int src, int dest, int cost) {
        adj[src].push_back(make_pair(dest, cost));
        cout << "new edge from:," << src << ",dest:," << dest << cost << endl;
    }

    void findPaths(int start, int goal, int budget) {
        vector<int> path;
        int pathCost = 0;
        ofstream outFile("paths.txt"); // Open the output file here to pass it to the utility function
        findPathsUtil(start, goal, budget, path, pathCost, outFile);
        outFile.close(); // Close the file after all paths are written
    }

    void loadEdgesFromFile(const string& filename) {
        ifstream file(filename);
        string line;
        while (getline(file, line)) {
            stringstream ss(line);
            int src, dest, cost;
            char delim;
            ss >> src >> delim >> dest >> delim >> cost;
            addEdge(src, dest, cost);
        }
    }

private:
    void findPathsUtil(int currNode, int goal, int budget, vector<int>& path, int& pathCost, ofstream& outFile) {
        path.push_back(currNode);
        cout << "currNode:," << currNode << endl;

        if (currNode == goal && pathCost <= budget) {
            cout << "goal found" << endl;
            printPath(path, pathCost, outFile);
        }

        for (auto i = adj[currNode].begin(); i != adj[currNode].end(); ++i) {
            int nextNode = i->first;
            int nextCost = i->second;
            cout << "nextNode:" << i->first << ",nextCost:," << i->second << endl;

            if (pathCost + nextCost <= budget) {
                pathCost += nextCost;
                findPathsUtil(nextNode, goal, budget, path, pathCost, outFile);
                pathCost -= nextCost; // Backtrack
            }
        }

        path.pop_back();
    }

    void printPath(vector<int>& path, int cost, ofstream& outFile) {
        for (int i = 0; i < path.size(); ++i) {
            outFile << path[i];
            if (i != path.size() - 1) outFile << ",";
        }
        //outFile << ", | Cost: " << cost << endl; // Include the cost at the end of each line
        outFile << endl;
    }
};

int main() {
    Graph g;
    
    string filename;
    cout << "Enter CSV filename: ";
    cin >> filename;

    g.loadEdgesFromFile(filename);

    int budget, start, goal;
    cout << "Enter budget, start node, goal node: ";
    cin >> budget >> start >> goal;
    
    cout<< "budget:," <<budget<<",start:,"<< start<< ",goal:,"<< goal<<endl;

    g.findPaths(start, goal, budget);

    cout << "Paths have been written to paths.txt" << endl;

    return 0;
}
