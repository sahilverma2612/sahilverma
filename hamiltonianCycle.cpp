#include <iostream>
#include <vector>
using namespace std;

#define V 5  // Change this based on number of vertices

// Utility function to check if current vertex can be added to the Hamiltonian Cycle
bool isSafe(int v, vector<vector<int>>& graph, vector<int>& path, int pos) {
    // Check if current vertex is adjacent to the previous vertex
    if (graph[path[pos - 1]][v] == 0)
        return false;

    // Check if the vertex has already been included
    for (int i = 0; i < pos; i++)
        if (path[i] == v)
            return false;

    return true;
}

// Recursive utility to solve the Hamiltonian Cycle problem
bool hamiltonianCycleUtil(vector<vector<int>>& graph, vector<int>& path, int pos) {
    // Base case: all vertices are in the path
    if (pos == V) {
        // Check if there is an edge from the last vertex to the first
        return graph[path[pos - 1]][path[0]] == 1;
    }

    // Try different vertices
    for (int v = 1; v < V; v++) {
        if (isSafe(v, graph, path, pos)) {
            path[pos] = v;

            if (hamiltonianCycleUtil(graph, path, pos + 1))
                return true;

            // Backtrack
            path[pos] = -1;
        }
    }
    return false;
}

bool hamiltonianCycle(vector<vector<int>>& graph) {
    vector<int> path(V, -1);
    path[0] = 0;  // Start from vertex 0

    if (!hamiltonianCycleUtil(graph, path, 1)) {
        cout << "No Hamiltonian Cycle exists\n";
        return false;
    }

    cout << "Hamiltonian Cycle found:\n";
    for (int i = 0; i < V; i++)
        cout << path[i] << " ";
    cout << path[0] << endl;  // Return to the start
    return true;
}

int main() {
    vector<vector<int>> graph = {
        {0, 1, 0, 1, 0},
        {1, 0, 1, 1, 1},
        {0, 1, 0, 0, 1},
        {1, 1, 0, 0, 1},
        {0, 1, 1, 1, 0}
    };

    hamiltonianCycle(graph);
    return 0;
}
