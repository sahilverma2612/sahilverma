#include <iostream>
#include <vector>
#include <string>
using namespace std;

// Function to compute LCS length
int lcs(string X, string Y) {
    int m = X.length();
    int n = Y.length();

    // Create a 2D DP table
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    // Build the dp table from bottom-up
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (X[i - 1] == Y[j - 1])
                dp[i][j] = 1 + dp[i - 1][j - 1];
            else
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
        }
    }

    return dp[m][n];
}

// Function to reconstruct and print the LCS string (optional)
string getLCS(string X, string Y) {
    int m = X.length(), n = Y.length();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    // Fill dp table
    for (int i = 1; i <= m; i++)
        for (int j = 1; j <= n; j++)
            if (X[i - 1] == Y[j - 1])
                dp[i][j] = 1 + dp[i - 1][j - 1];
            else
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);

    // Reconstruct LCS string from dp table
    string lcsStr;
    int i = m, j = n;
    while (i > 0 && j > 0) {
        if (X[i - 1] == Y[j - 1]) {
            lcsStr = X[i - 1] + lcsStr;
            i--;
            j--;
        } else if (dp[i - 1][j] > dp[i][j - 1])
            i--;
        else
            j--;
    }

    return lcsStr;
}

int main() {
    string str1, str2;
    cout << "Enter first string: ";
    cin >> str1;

    cout << "Enter second string: ";
    cin >> str2;

    int length = lcs(str1, str2);
    string sequence = getLCS(str1, str2);

    cout << "Length of LCS: " << length << endl;
    cout << "LCS: " << sequence << endl;

    return 0;
}
