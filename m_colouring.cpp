#include <iostream>
#include <vector>
using namespace std;

bool isSafe(int node, int color[], vector<vector<int>>& graph, int n, int col) {
    for (int k = 0; k < n; k++) {
        if (graph[node][k] == 1 && color[k] == col)
            return false;
    }
    return true;
}

bool solve(int node, int color[], int m, int n, vector<vector<int>>& graph) {
    if (node == n) return true;

    for (int col = 1; col <= m; col++) {
        if (isSafe(node, color, graph, n, col)) {
            color[node] = col;
            if (solve(node + 1, color, m, n, graph)) return true;
            color[node] = 0; // backtrack
        }
    }

    return false;
}

void graphColoring(vector<vector<int>>& graph, int m) {
    int n = graph.size();
    int color[n];
    for (int i = 0; i < n; i++) color[i] = 0;

    if (solve(0, color, m, n, graph)) {
        cout << "Solution exists. Vertex colors are:\n";
        for (int i = 0; i < n; i++) {
            cout << "Vertex " << i << " ---> Color " << color[i] << endl;
        }
    } else {
        cout << "No solution exists with " << m << " colors.\n";
    }
}

int main() {
    int n, m;
    cout << "Enter number of vertices: ";
    cin >> n;

    vector<vector<int>> graph(n, vector<int>(n, 0));
    cout << "Enter adjacency matrix:\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> graph[i][j];

    cout << "Enter number of colors (M): ";
    cin >> m;

    graphColoring(graph, m);

    return 0;
}
