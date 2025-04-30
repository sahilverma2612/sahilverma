#include <iostream>
#include <vector>
#include <climits>
using namespace std;

int tsp(vector<vector<int>>& graph, vector<bool>& visited, int pos, int n, int count, int cost, int start) {
    if (count == n && graph[pos][start]) {
        return cost + graph[pos][start]; // Return to the start
    }

    int ans = INT_MAX;

    for (int i = 0; i < n; i++) {
        if (!visited[i] && graph[pos][i]) {
            visited[i] = true;
            int temp = tsp(graph, visited, i, n, count + 1, cost + graph[pos][i], start);
            ans = min(ans, temp);
            visited[i] = false; // Backtrack
        }
    }

    return ans;
}

int main() {
    int n;
    cout << "Enter number of cities: ";
    cin >> n;

    vector<vector<int>> graph(n, vector<int>(n));
    cout << "Enter the cost adjacency matrix:\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> graph[i][j];

    vector<bool> visited(n, false);
    visited[0] = true; // Starting from city 0

    int result = tsp(graph, visited, 0, n, 1, 0, 0);

    cout << "Minimum travelling cost: " << result << endl;
    return 0;
}
