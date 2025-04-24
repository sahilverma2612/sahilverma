**1. Selection Sorting:-** 

Time Complexity Analysis:-
This C++ program measures the average execution time of the Selection Sort algorithm 
for arrays of increasing size, ranging from 1,000 to 25,000 elements. 
It helps analyze how the time complexity of the algorithm grows with the input size. Inefficient for large datasets.

  Time Complexity:

    Best Case: O(n²)

    Average Case: O(n²)

    Worst Case: O(n²)

#include <bits/stdc++.h>

using namespace std;

void selectionSort(long int arr[], int size) {

    for (int i = 0; i < size - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < size; j++) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        if (minIndex != i) {
            swap(arr[i], arr[minIndex]);
        }
    }
}

int main() {

    srand(time(0));
    int n = 1000;
    while (n <= 25000) {
        clock_t totalTime = 0;
        for (int i = 1; i <= 10; i++) {
            long int arr[n];
            for (int j = 0; j < n; j++) {
                arr[j] = rand() % n + 1;
            }
            clock_t start = clock();
            selectionSort(arr, n);
            clock_t end = clock();
            totalTime += (end - start);
        }
        double avgTime = (double)totalTime / (10 * CLOCKS_PER_SEC);
        cout << "Size: " << n << ", Time: " << avgTime << " seconds" << endl;
        n += 1000;
    }
    return 0;
}

**2. Bubble Sort:-**

Time Complexity Analysis:-This C++ program measures the average execution time of the Bubble Sort algorithm 
for arrays of increasing size, from 1,000 to 25,000 elements. It demonstrates how the time complexity 
of Bubble Sort grows with larger input sizes.
**Time Complexity:**

    Best Case: O(n) (only when the array is already sorted)

    Average Case: O(n²)

    Worst Case: O(n²)

#include <bits/stdc++.h>

  using namespace std;

void bubbleSort(long int arr[], int n) {

    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1 - i; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(arr[j], arr[j + 1]);
            }
        }
    }
}

int main() {

    srand(time(0));
    int n = 1000;
    while (n <= 25000) {
        clock_t totalTime = 0;
        for (int i = 1; i <= 10; i++) {
            long int arr[n];
            for (int j = 0; j < n; j++) {
                arr[j] = rand() % n + 1;
            }
            clock_t start = clock();
            bubbleSort(arr, n);
            clock_t end = clock();
            totalTime += (end - start);
        }
        double avgTime = (double)totalTime / (10 * CLOCKS_PER_SEC);
        cout << "Size: " << n << ", Time: " << avgTime << " seconds" << endl;
        n += 1000;
    }
    return 0;
}

**3. Quick Sort:-**

    Time Complexity Analysis:-This C++ program measures the average execution time of the Quick Sort algorithm 
    on arrays of increasing size, ranging from 1,000 to 25,000 elements. It helps visualize the efficiency
    and performance of Quick Sort in practical scenarios.

  **Case	Complexity:-**
 
        Best	           O(n log n)
        Average	           O(n log n)
        Worst	           O(n²)

#include <bits/stdc++.h>
using namespace std;

int partition(long int arr[], int low, int high) {

    long int pivot = arr[low];
     long int i = low;
    long int j = high;

    while (i < j) {
        while (arr[i] <= pivot && i <= high - 1) {
            i++;
        }

        while (arr[j] > pivot && j >= low + 1) {
            j--;
        }
        if (i < j) swap(arr[i], arr[j]);
    }
    swap(arr[low], arr[j]);
    return j;
}

void quickSort(long int arr[], int low, int high) {

    if (low < high) {
        int pIndex = partition(arr, low, high);
        quickSort(arr, low, pIndex - 1);
        quickSort(arr, pIndex + 1, high);
        
    }
}

int main() {

    srand(time(0)); // Seed for random numbers

    int n = 1000;
    while (n <= 25000) {
        double totalTime = 0.0;

        for (int i = 1; i <= 10; i++) {
            long int arr[n];
            
            // Generate random array of size n
            for (int j = 0; j < n; j++) {
                arr[j] = rand() % n + 1;
            }

            auto start = chrono::high_resolution_clock::now();
            quickSort(arr, 0, n - 1);
            auto end = chrono::high_resolution_clock::now();

            chrono::duration<double> duration = end - start;
            totalTime += duration.count();
        }

        double avgTime = totalTime / 10;
        cout << "Size: " << n << ", Time: " << avgTime << " seconds" << endl;

        n += 1000; // Increment n properly
    }

    return 0;
}   

**4.Merge Sort:-** 

Time Complexity Analysis:-This C++ program measures the average execution time of the Merge Sort algorithm 
for increasing input sizes. It demonstrates the efficiency and scalability of Merge Sort,
which is known for its reliable performance across diverse datasets.

Case	Complexity:-

    Best	        O(n log n)
    Average	        O(n log n)
    Worst	        O(n log n)

 #include <bits/stdc++.h>
 
  #include <ctime>
  using namespace std;

void merge(long int arr[], long int low, long int mid, long int high) {

    long int i = low;
    long int j = mid + 1;
    long int k = 0;
    vector<long int> temp(high - low + 1);

    while (i <= mid && j <= high) {
        if (arr[i] <= arr[j]) {
            temp[k] = arr[i];
            i++;
        } else {
            temp[k] = arr[j];
            j++;
        }
        k++;
    }

    while (i <= mid) {
        temp[k++] = arr[i++];
    }

    while (j <= high) {
        temp[k++] = arr[j++];
    }

    for (long int i = low, k = 0; i <= high; i++, k++) {
        arr[i] = temp[k];
    }
}

void mergeSort(long int arr[], long int low, long int high) {

    if (low >= high) return;
    long int mid = (low + high) / 2;
    mergeSort(arr, low, mid);
    mergeSort(arr, mid + 1, high);
    merge(arr, low, mid, high);
}

int main() {

    long int n = 10000;
    int it = 0;
    double time[10];

    cout << "size, time\n";

    while (it < 10) {
        long int a[n];
        for (int i = 0; i < n; i++) {
            a[i] = rand() % n + 1;
        }

        clock_t start = clock();
        mergeSort(a, 0, n - 1);
        clock_t end = clock();

        time[it] = double(end - start) / CLOCKS_PER_SEC;

        cout << n << ", " << time[it] << "," << endl;

        n += 10000;
        it++;
    }

    return 0;
}

**5. Iterative QuickSort (QSI):-**

    This program implements the Iterative QuickSort algorithm to sort an array of integers
    using the divide and conquer technique. The code eliminates the need for recursion
    by utilizing a stack, making it suitable for large arrays where recursion might
    cause stack overflow.

  **Time Complexity:-**
  
    Average Case: O(n log n)

    Worst Case: O(n²)

    #include<bits/stdc++.h>
    using namespace std;

    int Partition(long int arr[], long int low, long int high)
    {
    long int pivot = arr[low];
    long int i = low + 1;
    long int j = high;
    while (i <= j)
    {
        while (i <= j && arr[i] <= pivot)
        {
            i++;
        }
        while (i <= j && arr[j] > pivot)
        {
            j--;
        }
        if (i < j)
        {
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[low], arr[j]);
    return j;
    }

    void QSI(long int arr[], int low, int high)
    {
    stack<int> st;
    do {
        while (low < high)
        {
            int j = Partition(arr, low, high);
            if ((j - low) > (high - j))
            {
                st.push(low);
                st.push(j - 1);
                low = j + 1;
            }
            else
            {
                st.push(j + 1);
                st.push(high);
                high = j - 1;
            }
        }
        if (st.empty()) return;
        high = st.top();
        st.pop();
        low = st.top();
        st.pop();
    } while (true);
    }

    bool isSorted(long int arr[], int n) {
    for (int i = 1; i < n; i++) {
        if (arr[i] < arr[i - 1]) return false;
    }
    return true;
    }

    int main()
    {
    srand(time(0));

    int n = 1000;
    while (n <= 10000) {
        double totalTime = 0.0;

        for (int i = 1; i <= 10; i++) {
            long int arr[n];
            
            for (int j = 0; j < n; j++) {
                arr[j] = rand() % n + 1;
            }

            clock_t start = clock();
            QSI(arr, 0, n - 1);
            clock_t end = clock();

            double duration = double(end - start) / CLOCKS_PER_SEC;
            totalTime += duration;

            if (!isSorted(arr, n)) {
                cout << "Sorting failed for size " << n << endl;
                return 1;
            }
        }

        double avgTime = totalTime / 10;
        cout << "Size: " << n << ", Time: " << avgTime << " seconds" << endl;

        n += 1000;
    }

    return 0;
    }
**6.Insertion Sorting:-**

      It is a simple sorting algorithm that builds the final sorted array one element at a time. 
      It is much like sorting playing cards in your hands—picking one card at a time 
      and placing it in the correct position.

   **Time Complexity:-**

      | Case        | Complexity |
      |-------------|------------|
      | Best Case   | O(n)       |
      | Average Case| O(n²)      |
      | Worst Case  | O(n²)      |

    #include <iostream>
    #include <vector>
    #include <cstdlib>
    #include <ctime>

    using namespace std;

    int main() 
    { 
    srand(time(0));
    int n = 1000; 

    while (n <= 10000) 
    { 
        long long total_time = 0; 

        for (int i = 1; i <= 10; i++) 
        { 
            vector<int> arr(n); 
            for (int j = 0; j < n; j++) 
            { 
                arr[j] = rand(); 
            } 

            clock_t start_time = clock(); 

            for (int i = 1; i < n; i++) 
            { 
                int j = i - 1; 
                int temp = arr[i]; 
                while (j >= 0 && arr[j] > temp) 
                { 
                    arr[j + 1] = arr[j]; 
                    j = j - 1; 
                } 
                arr[j + 1] = temp; 
            } 

            clock_t end_time = clock(); 
            long long duration = (end_time - start_time) * 1000000 / CLOCKS_PER_SEC; 

            total_time += duration; 
        } 

        double avg_time = total_time / 10.0 / 1000.0; 

        cout << "Size of array: " << n << " | Avg Time (ms): " << avg_time << endl; 
        
        n += 1000; 
    } 

    return 0;  
    }

    # 🧳 Fractional Knapsack Problem in C++

**7. Knapsack Problem:-**

This project implements the **Fractional Knapsack Problem** using a greedy approach in C++. In this variant of the classic
knapsack problem, we are allowed to take **fractions of items** to maximize the total profit, given a fixed weight capacity.


## 🚀 Features

- Sorts items based on **profit-to-weight ratio**.
- Selects the best items (or their fractions) to maximize profit.
- Handles **full and partial item selection**.
- Simple and efficient greedy algorithm.

---

## 💡 Problem Statement

Given a list of items with profits and weights, and a knapsack with a limited capacity, maximize
the total profit by taking whole or fractional parts of items.



## 📈 Time Complexity

- Sorting: **O(n log n)** (based on profit/weight ratio)
- Selection: **O(n)**
- **Total:** O(n log n)

      #include<bits/stdc++.h>
      using namespace std;

      struct Item{
      int profit,weight;
};

    bool compare(Item a, Item b ){
    double r1 = (double)a.profit/a.weight;
    double r2 = (double)b.profit/b.weight;
    return r1>r2;
    }

    void knapsack(vector<Item>Items,int capacity){
    //sort kar do profit/weight ratio ke accoding decreasing order me 
    sort(Items.begin(),Items.end(),compare);

    double totalprofit = 0.0;
    int currentweight =0;
     for(int i =0;i<Items.size();i++){
      if(currentweight + Items[i].weight <= capacity ){
        //pura item le lo
        totalprofit = totalprofit + Items[i].profit;
        currentweight = currentweight + Items[i].weight;
        cout << "Profit: " << Items[i].profit << ", Weight: " << Items[i].weight << " (Full)\n";
      }
      else{
        //fractional part le lo agar pura nahi le sakte to
        int remainweight = capacity - currentweight;
        double fraction_profit = Items[i].profit * ((double)remainweight / Items[i].weight);
        totalprofit = totalprofit + fraction_profit;
        cout << "Profit: " << fraction_profit << ", Weight: " << remainweight << " (Fraction)\n";
        break;//knapsack full ho gya
      }
     }
     cout << "Total Profit: " << totalprofit << endl;
    }

    int main(){
    vector<Item> Items = {{25,18},{24,15},{15,10}};

    knapsack(Items,20);
}

**8. Magic Square:-**


This program generates a **Magic Square** of any given odd order `n` using the **Siamese method** (also known as the De la Loubère's method).
A magic square is a grid of distinct positive integers such that the sum of the integers in each row, column, and both main diagonals 
is the same.


## ✨ Features

- Generates a **Magic Square for any odd number** `n`.
- Uses a simple and efficient algorithm.
- Output is formatted in a readable grid.


## 🧠 How It Works (Siamese Method)

1. Start with 1 in the middle of the top row.
2. Place the next number one row up and one column right.
3. If that cell is already filled or out of bounds:
   - Wrap around the grid using modulo.
   - If still filled, move down one cell instead.
4. Repeat until the square is full.

---

## 📈 Time Complexity

- **O(n²)** for filling the entire grid.

      #include <bits/stdc++.h>
       using namespace std;

      void magicsquare(vector<vector<int>> &a, int n) {
      int i = 0, j = n / 2; // Start at the middle of the first row

      for (int k = 1; k <= n * n; k++) {
        a[i][j] = k; // Place the number in the grid

        // Store current position before moving
        i = (i - 1 + n) % n; // Move up (wrap around using modulo)
        j = (j + 1) % n;     // Move right (wrap around using modulo)

        // If the new cell is already occupied, move down instead
        if (a[i][j] != 0) {
            i = (i + 1) % n; // Move down
        
        }
      } 
    }

      int main() {
      int n;
      cin >> n;

      if (n % 2 == 0) {
        cout << "Magic square is only possible for odd numbers." << endl;
        return 1;
      }

      // Creating a 2D vector initialized with 0
      vector<vector<int>> a(n, vector<int>(n, 0));

      // Generate magic square
      magicsquare(a, n);

      // Display the magic square
      cout << "Magic Square of order " << n << ":\n";
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(4) << a[i][j]; // Print formatted numbers
        }
        cout << endl;
      }

      return 0;
      }

 **Activity Selection Problem**
 
This program solves the Activity Selection Problem using a greedy algorithm. The goal is to
 select the maximum number of activities that can be scheduled in a given time period
such that no two activities overlap.

**Problem Description**
Given a set of activities, each having a start time and a finish time, the task is to 
select the maximum number of activities that don't overlap with each other. 
The activities must be selected such that for any two selected activities,
the start time of one activity should be greater than or equal to the
finish time of the previous selected activity.

    #include<bits/stdc++.h>
    using namespace std;

    struct Activity {
    int start, finish;
    };

    bool compare(Activity a1, Activity a2) {
    return a1.finish < a2.finish;
}

    void activitySelection(vector<Activity> activities) {
    sort(activities.begin(), activities.end(), compare);

    cout << "Selected Activities: \n";
    int lastFinishTime = activities[0].finish;
    cout << "(" << activities[0].start << ", " << activities[0].finish << ")\n";
    
    for (int i = 1; i < activities.size(); i++) {
        if (activities[i].start >= lastFinishTime) {
            cout << "(" << activities[i].start << ", " << activities[i].finish << ")\n";
            lastFinishTime = activities[i].finish;
        }
    }
}

    int main() {
    vector<Activity> activities = {{1, 3}, {2, 5}, {3, 9}, {0, 6}, {5, 7}, {8, 9}};
    
    activitySelection(activities);
    
    return 0;
}



**Peak Element in an Array**

This program finds the peak element in an array using an optimized binary search approach.

**Problem Description**

A peak element in an array is an element that is greater than or equal to its neighbors. 
For example, for an array [1, 3, 20, 4, 1], 20 is a peak element as it is greater
than both its neighbors, 3 and 4.


**Time Complexity**: O(n) [brute force solution]

Optimized Binary Search Approach (Used in the code):

This method uses binary search to efficiently find a peak element in O(log n) time.

Time Complexity: O(log n)


    #include<bits/stdc++.h>
    using namespace std;

    // Optimized approach using Binary search to find a peak element
    int peak(int* arr, int n) {
    if (n == 1) return 0;
    if (arr[0] > arr[1]) return 0;
    if (arr[n - 1] > arr[n - 2]) return n - 1;
    
    int low = 1, high = n - 2;
    while (low <= high) {
        int mid = (low + high) / 2;
        if ((arr[mid] > arr[mid - 1]) && (arr[mid] > arr[mid + 1])) return mid;
        else if (arr[mid] > arr[mid - 1]) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return -1;
    }

    int main() {
    int arr[5] = {11, 13, 20, 25, 18};
    int res = peak(arr, 5);
    cout << "Peak element is: " << arr[res] << endl;
    return 0;
    }



**2D Peak Element**


This program finds a peak element in a 2D matrix using a binary search-based approach.
A peak element is an element that is strictly greater than its neighbors in the matrix.


    #include<bits/stdc++.h>
    using namespace std;

    int findmax(vector<vector<int>>& mat, int n, int m, int col) {
    int max = -1;
    int max_index = -1;
    for (int i = 0; i < n; i++) {
        if (mat[i][col] > max) {
            max = mat[i][col];
            max_index = i;
        }
    }
    return max_index;
    }

    vector<int> findpeak(vector<vector<int>>& mat) {
    int n = mat.size();
    int m = mat[0].size();
    int low = 0;
    int high = m - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        int max_index = findmax(mat, n, m, mid);
        if ((mid == 0 || mat[max_index][mid - 1] < mat[max_index][mid]) && (mid == m - 1 
        || mat[max_index][mid] > mat[max_index][mid + 1])) {
            return {max_index, mid};
        } else if (mid > 0 && mat[max_index][mid - 1] > mat[max_index][mid]) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }
    return {-1, -1};
    }

    int main() {
    vector<vector<int>> mat = {
        {10, 8, 10, 10},
        {14, 13, 12, 11},
        {15, 9, 11, 21},
        {16, 17, 19, 20}
    };
    vector<int> res = findpeak(mat);
    cout << "Peak element is: " << mat[res[0]][res[1]] << endl;
    return 0;
    }

    vector<vector<int>> mat = {
    {10, 8, 10, 10},
    {14, 13, 12, 11},
    {15, 9, 11, 21},
    {16, 17, 19, 20}
    };

  **Output:-**
        Peak element is: 21


**Prim’s Algorithm Using Near Array for Minimum Spanning Tree (MST)**


This program implements Prim’s Algorithm using the Near Array technique to find 
the Minimum Spanning Tree (MST) of a graph. The graph is represented using an 
adjacency matrix where the element at cost[i][j] represents the weight of the 
edge between vertex i and vertex j.


**Problem Description**

Given a graph with n vertices and an adjacency matrix representing the edge weights,
the task is to find the Minimum Spanning Tree (MST). The MST is a subset of the 
edges that connect all the vertices without any cycles and with the minimum
possible total edge weight.


    #include <iostream>
    #include <vector>
    #include <climits>

    using namespace std;

    void primsNear(int n, vector<vector<int>> &cost) {
    vector<int> near(n, INT_MAX);  // Near array
    vector<pair<int, int>> MST;    // Store MST edges
    int minCost = 0;

    // Step 1: Find the first minimum cost edge
    int u = -1, v = -1, minEdge = INT_MAX;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (cost[i][j] < minEdge) {
                minEdge = cost[i][j];
                u = i;
                v = j;
            }
        }
    }

    MST.push_back({u, v});
    minCost += minEdge;
    
    // Step 2: Initialize Near Array
    for (int i = 0; i < n; i++) {
        if (cost[i][u] < cost[i][v]) {
            near[i] = u;
        } else {
            near[i] = v;
        }
    }
    near[u] = near[v] = -1; // Mark included vertices

    // Step 3: Add remaining n-2 edges
    for (int i = 1; i < n - 1; i++) {
        int next = -1, minVal = INT_MAX;

        // Find the next minimum cost edge
        for (int j = 0; j < n; j++) {
            if (near[j] != -1 && cost[j][near[j]] < minVal) {
                minVal = cost[j][near[j]];
                next = j;
            }
        }

        MST.push_back({next, near[next]});
        minCost += cost[next][near[next]];
        near[next] = -1; // Mark included vertex

        // Update Near Array
        for (int j = 0; j < n; j++) {
            if (near[j] != -1 && cost[j][next] < cost[j][near[j]]) {
                near[j] = next;
            }
        }
    }

    // Print the MST
    cout << "Minimum Spanning Tree (MST):\n";
    for (int i = 0; i < MST.size(); i++) {
        cout << MST[i].first << " - " << MST[i].second << "\n";
    }
    cout << "Total Minimum Cost: " << minCost << endl;
    }

    int main() {
    int n = 5; // Number of vertices
    vector<vector<int>> cost = {
        {INT_MAX, 2, INT_MAX, 6, INT_MAX},
        {2, INT_MAX, 3, 8, 5},
        {INT_MAX, 3, INT_MAX, INT_MAX, 7},
        {6, 8, INT_MAX, INT_MAX, 9},
        {INT_MAX, 5, 7, 9, INT_MAX}
    };

    primsNear(n, cost);
    return 0;
    }

INPUT:-

    vector<vector<int>> cost = {
    {INT_MAX, 2, INT_MAX, 6, INT_MAX},
    {2, INT_MAX, 3, 8, 5},
    {INT_MAX, 3, INT_MAX, INT_MAX, 7},
    {6, 8, INT_MAX, INT_MAX, 9},
    {INT_MAX, 5, 7, 9, INT_MAX}
    };

Output:-

        Minimum Spanning Tree (MST):
        0 - 1
        1 - 2
        1 - 4
        0 - 3
        Total Minimum Cost: 16



   **Strassen's Matrix Multiplication Algorithm**
   
This program implements Strassen's Algorithm for matrix multiplication, which is an optimized algorithm
to multiply two square matrices in a time complexity of 𝑂(𝑛log27) compared to the
conventional 𝑂(𝑛3)  matrix multiplication.The algorithm divides the input 
matrices into smaller sub-matrices and performs matrix multiplications recursively.

**Steps Followed in Strassen's Algorithm:**

1.  Divide the matrices into smaller sub-matrices (quadrants).

2.  Compute 7 intermediate products (P1 to P7) using a combination of matrix additions and subtractions.

3.  Combine the sub-matrices to get the resulting matrix.



         #include <iostream>
        #include <vector>

        using namespace std;

        typedef vector<vector<int>> Matrix;

        // Function to add two matrices
        Matrix addMatrix(const Matrix &A, const Matrix &B, int size) {
        Matrix result(size, vector<int>(size, 0));
        for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            result[i][j] = A[i][j] + B[i][j];
        return result;
        }

        // Function to subtract two matrices
        Matrix subtractMatrix(const Matrix &A, const Matrix &B, int size) {
         Matrix result(size, vector<int>(size, 0));
        for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            result[i][j] = A[i][j] - B[i][j];
         return result;
         }

        // Strassen’s Algorithm for matrix multiplication
        Matrix strassenMultiply(Matrix A, Matrix B, int size) {
        if (size == 1) {
        return {{A[0][0] * B[0][0]}};
        }

        int newSize = size / 2;
        Matrix A11(newSize, vector<int>(newSize)), A12(newSize, vector<int>(newSize));
        Matrix A21(newSize, vector<int>(newSize)), A22(newSize, vector<int>(newSize));
        Matrix B11(newSize, vector<int>(newSize)), B12(newSize, vector<int>(newSize));
        Matrix B21(newSize, vector<int>(newSize)), B22(newSize, vector<int>(newSize));

        // Dividing matrices into quadrants
        for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + newSize];
            A21[i][j] = A[i + newSize][j];
            A22[i][j] = A[i + newSize][j + newSize];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + newSize];
            B21[i][j] = B[i + newSize][j];
            B22[i][j] = B[i + newSize][j + newSize];
        }
        }

        // Computing Strassen's 7 multiplications
        Matrix P1 = strassenMultiply(A11, subtractMatrix(B12, B22, newSize), newSize);
        Matrix P2 = strassenMultiply(addMatrix(A11, A12, newSize), B22, newSize);
        Matrix P3 = strassenMultiply(addMatrix(A21, A22, newSize), B11, newSize);
        Matrix P4 = strassenMultiply(A22, subtractMatrix(B21, B11, newSize), newSize);
        Matrix P5 = strassenMultiply(addMatrix(A11, A22, newSize), addMatrix(B11, B22, newSize), newSize);
        Matrix P6 = strassenMultiply(subtractMatrix(A12, A22, newSize), addMatrix(B21, B22, newSize), newSize);
        Matrix P7 = strassenMultiply(subtractMatrix(A11, A21, newSize), addMatrix(B11, B12, newSize), newSize);

        // Computing the resultant quadrants
        Matrix C11 = addMatrix(subtractMatrix(addMatrix(P5, P4, newSize), P2, newSize), P6, newSize);
        Matrix C12 = addMatrix(P1, P2, newSize);
        Matrix C21 = addMatrix(P3, P4, newSize);
        Matrix C22 = subtractMatrix(subtractMatrix(addMatrix(P5, P1, newSize), P3, newSize), P7, newSize);

        // Combining the quadrants into a single matrix
        Matrix C(size, vector<int>(size, 0));
        for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            C[i][j] = C11[i][j];
            C[i][j + newSize] = C12[i][j];
            C[i + newSize][j] = C21[i][j];
            C[i + newSize][j + newSize] = C22[i][j];
        }
        }

        return C;
        }

        // Function to print a matrix
        void printMatrix(const Matrix &M) {
        for (const auto &row : M) {
        for (int val : row)
            cout << val << " ";
        cout << endl;
        }
        }

        int main() {
        int n = 4;  // Matrix size (should be a power of 2)
        Matrix A = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        Matrix B = {{17, 18, 19, 20}, {21, 22, 23, 24}, {25, 26, 27, 28}, {29, 30, 31, 32}};

        cout  << "Matrix A:\n";
        printMatrix(A);
        cout << "\nMatrix B:\n";
        printMatrix(B);

        Matrix C = strassenMultiply(A, B, n);

        cout << "\nResultant Matrix (A * B):\n";
        printMatrix(C);

        return 0;
          }

    
    INPUT:-
        
        Matrix A = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        Matrix B = {{17, 18, 19, 20}, {21, 22, 23, 24}, {25, 26, 27, 28}, {29, 30, 31, 32}};

    Output:-

        Matrix A:
        1 2 3 4 
        5 6 7 8 
        9 10 11 12 
        13 14 15 16 

        Matrix B:
        17 18 19 20 
        21 22 23 24 
        25 26 27 28 
        29 30 31 32 

        Resultant Matrix (A * B):
        250 260 270 280 
        570 596 622 648 
        890 932 974 1016 
        1210 1268 1326 1384 

**Kruskal’s Algorithm for Minimum Spanning Tree (MST)**

This program implements Kruskal's Algorithm for finding the Minimum Spanning Tree (MST)
of a weighted undirected graph. The algorithm uses the Union-Find (Disjoint Set Union) 
data structure with path compression and union by size for efficient cycle detection 
and edge merging.


    #include<bits/stdc++.h>
    using namespace std;

    class DSU{
    vector<int>p;
    public:
    DSU(int n)
    {   p.resize(n);
        for(int i=0;i<n;i++)
        {
            p[i]=-1;
        }
    }

    int find(int u) {
        if (p[u] < 0) return u;  
        return p[u] = find(p[u]);  // Path compression
    }
    
    //Union by size
    bool unionSet(int u, int v) {
        int rootU = find(u);
        int rootV = find(v);
                
        if (rootU == rootV) return false; // Already connected

        if (p[rootU] < p[rootV]) {  // rootU ka size bada hai
            p[rootU] += p[rootV];
            p[rootV] = rootU;
        } else {
            p[rootV] += p[rootU];
            p[rootU] = rootV;
        }
        return true;
    }
    };

    // Custom Comparator for Min-Heap
    struct Compare {
    bool operator()(vector<int>& a, vector<int>& b) {
        return a[0] > b[0]; // Min-Heap based on weight
    }
    };

      // Kruskal’s Algorithm for MST
      void kruskalMST(int n, vector<vector<int>>& edges) {
    priority_queue<vector<int>, vector<vector<int>>, Compare> minHeap(edges.begin(), edges.end());

    DSU dsu(n);
    int minCost = 0, edgesUsed = 0;

    while(!minHeap.empty() && edgesUsed < n - 1) {
        auto edge = minHeap.top();
        minHeap.pop();

        int w = edge[0], u = edge[1], v = edge[2];

        if(dsu.unionSet(u, v)) { // If adding this edge doesn’t create a cycle
            minCost += w;
            edgesUsed++;
            cout << "Edges: (" << u << "-" << v << ") with cost: " << w << endl;
        }
    }

    if (edgesUsed == n - 1)
        cout << "Minimum Spanning Tree Cost: " << minCost << endl;
    else
        cout << "MST not possible!" << endl;
    }

    int main() {
    int n = 5; // Number of vertices
    vector<vector<int>> edges = {
        {2, 0, 1}, {3, 1, 2}, {7, 2, 4}, {5, 1, 4}, {9, 3, 4}, {6, 0, 3}, {8, 1, 3}
    };

    kruskalMST(n, edges);

    return 0;
    }

  **INPUT:-**

    int n = 5; // Number of vertices
    vector<vector<int>> edges = {
    {2, 0, 1}, {3, 1, 2}, {7, 2, 4}, {5, 1, 4}, {9, 3, 4}, {6, 0, 3}, {8, 1, 3}
    };

**OUTPUT:-**

    Edges: (0-1) with cost: 2
    Edges: (1-2) with cost: 3
    Edges: (0-3) with cost: 6
    Edges: (1-4) with cost: 5
    Minimum Spanning Tree Cost: 16


  **Dijkstra's Algorithm – Shortest Path Finder**
  
This C++ program implements Dijkstra’s Algorithm to find the shortest path from a source vertex
to all other vertices in a weighted graph represented by an adjacency matrix.

🧠 Overview
Dijkstra's algorithm is a greedy algorithm used to solve the single-source shortest path problem
in graphs with non-negative edge weights. It efficiently computes the shortest distance from a 
source node to all other nodes in the graph.


    #include<bits/stdc++.h>
    using namespace std;

    #define V 9

    int minDistance(int dist[], bool sptSet[])
    {
	
	int min = INT_MAX, min_index;

	for (int v = 0; v < V; v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
    }


    void printSolution(int dist[], int n)
    {
	printf("Vertex Shortest Distance from Source\n");
	for (int i = 0; i < V; i++)
		printf("\t%d \t %d\n", i, dist[i]);
    }


    void dijkstra(int graph[V][V], int src)
    {
	int dist[V]; 

	bool sptSet[V]; 
	for (int i = 0; i < V; i++)
		dist[i] = INT_MAX, sptSet[i] = false;

	
	dist[src] = 0;

	
	for (int count = 0; count < V - 1; count++) {
		
 		int u = minDistance(dist, sptSet);

		
		sptSet[u] = true;

		
		for (int v = 0; v < V; v++)

			
			if (!sptSet[v] && graph[u][v]&& dist[u] != INT_MAX&& dist[u] + graph[u][v] < dist[v])
             dist[v] = dist[u] + graph[u][v];
	}

	
	printSolution(dist, V);
    }


    int main()
    {
	
	int graph[V][V] = { { 0, 4, 0, 0, 0, 0, 0, 8, 0 },
						{ 4, 0, 8, 0, 0, 0, 0, 11, 0 },
						{ 0, 8, 0, 7, 0, 4, 0, 0, 2 },
						{ 0, 0, 7, 0, 9, 14, 0, 0, 0 },
						{ 0, 0, 0, 9, 0, 10, 0, 0, 0 },
						{ 0, 0, 4, 14, 10, 0, 2, 0, 0 },
						{ 0, 0, 0, 0, 0, 2, 0, 1, 6 },
						{ 8, 11, 0, 0, 0, 0, 1, 0, 7 },
						{ 0, 0, 2, 0, 0, 0, 6, 7, 0 } };

	dijkstra(graph, 0);

	return 0;
    }

**EXAMPLE GRAPH:-**

    (0)--4--(1)--8--(2)--7--(3)
     |       |       |       |
     8       11      2       14
     |       |       |       |
    (7)--1--(6)--2--(5)--10--(4)

  **OUTPUT:-**
  
    Vertex Shortest Distance from Source
	0 	 0
	1 	 4
	2 	 12
	3 	 19
	4 	 21
	5 	 11
	6 	 9
	7 	 8
	8 	 14


 **Multi-Stage Graph**

 This C++ program finds the minimum-cost path from the source node to the destination 
 node in a directed acyclic multi-stage graph (DAG) using dynamic programming.

📘 Problem Overview
A multi-stage graph is a directed weighted graph where nodes are divided into different stages,
and all edges connect nodes from one stage to a later stage. The goal is to find the minimum-cost
path from the source (first stage) to the destination (last stage).


    #include<bits/stdc++.h>
    using namespace std;

    const int INF = 1e9;

    int countStages(vector<vector<int>>&Graph){
    int k = 1;//source node ka stage 1
    int u =0;
    for(int v = u+1;v<Graph.size();v++)
    {
    if(Graph[u][v]!=INF){
        u = v;
        k++;
    }
    }
    return k;
    }

    void MGP(vector<vector<int>>&Graph,int k){
    int n = Graph.size();
    vector<int> fdist(n, INF); // Minimum cost from each node to final node ko initialize kar do infinity se
    vector<int> d(n, -1);     // Stores next node in shortest path ko initialize kar do -1 sw

    fdist[n-1] = 0; //final node ka cost 0 hoga
    for(int j = n-2;j>=0;j--)
    {
        for(int r=j+1;r<n;r++)
        {
            if (Graph[j][r] != INF && Graph[j][r] != 0)// Valid edge
             { 
                if (Graph[j][r] + fdist[r] < fdist[j]) {
                    fdist[j] = Graph[j][r] + fdist[r];
                    d[j] = r; // next node ko store kar do  shortest path me 
                }
            }
        }
    }
   
    // Find minimum-cost path using d[]
    vector<int> p(k); // Path array
    p[0] = 1;  // First node (1-based indexing)
    p[k - 1] = n; // Last node (1-based indexing)

    for (int j = 1; j < k - 1; j++) {
    p[j] = d[p[j - 1] - 1] + 1; // Convert 0-based to 1-based
    }

    // Output the minimum-cost path
    cout << "\nMinimum-cost path: ";
    for (int i = 0; i < k; i++) {
    cout << p[i] << (i == k - 1 ? "\n" : " -> ");
    }


     cout << "Minimum cost: " << fdist[0] << endl;
    }
    int main()
    {   const int INF = 1e9;
    int n;
    cout << "Enter number of nodes: ";
    cin >> n;

    // Graph initialization: Diagonal = 0, Others = INF
    vector<vector<int>> Graph(n, vector<int>(n, INF));

    for (int i = 0; i < n; i++) {
        Graph[i][i] = 0; // Diagonal elements = 0
    }

    int edges;
    cout << "Enter number of edges: ";
    cin >> edges;

    cout << "Enter edges (u v cost) (1-based indexing):\n";
    for (int i = 0; i < edges; i++) {
        int u, v, cost;
        cin >> u >> v >> cost;
        Graph[u - 1][v - 1] = cost; // Convert to 0-based index
    }

    // Displaying the Graph
    cout << "\nGraph Representation (Adjacency Matrix):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (Graph[i][j] == INF)
                cout << "INF\t";
            else
                cout << Graph[i][j] << "\t";
        }
        cout << endl;
    }

     int k = countStages(Graph);
     cout<<"Number of Stages:"<<k<<endl;
     MGP(Graph,k);

     return 0;

    }
**INPUT:-**

    Enter number of nodes: 8
    Enter number of edges: 10
    Enter edges (u v cost) (1-based indexing):
    1 2 1
    1 3 2
    2 4 3
    2 5 4
    3 5 2
    4 6 6
    5 6 1
    5 7 3
    6 8 2
    7 8 1

**OUTPUT:-**

    Graph Representation (Adjacency Matrix):
    0	1	2	INF	INF	INF	INF	INF	
    INF	0	INF	3	4	INF	INF	INF	
    INF	INF	0	INF	2	INF	INF	INF	
    INF	INF	INF	0	INF	6	INF	INF	
    INF	INF	INF	INF	0	1	3	INF	
    INF	INF	INF	INF	INF	0	INF	2	
    INF	INF	INF	INF	INF	INF	0	1	
    INF	INF	INF	INF	INF	INF	INF	0	

    Number of Stages: 5

    Minimum-cost path: 1 -> 3 -> 5 -> 6 -> 8
    Minimum cost: 7


 **Floyd-Warshall Algorithm - All Pairs Shortest Path**
 
This C++ program implements the Floyd-Warshall Algorithm, which is used to 
find the shortest distances between every pair of vertices in a weighted
directed graph.


**Time Complexity: O(V³)
Space Complexity: O(V²)**


    #include<bits/stdc++.h> 
    using namespace std ; 
    void floydWarshall(vector<vector<int>> &adj, int n) 
    { 
    for(int k=0; k<n; k++) 
    { 
    for(int i=0; i<n; i++) 
    { 
    for(int j=0; j<n; j++) 
    { 
    if( ((adj[i][j] == -1) || (adj[i][j] > adj[i][k] + adj[k][j])) && (adj[k][j] 
    != -1 && adj[i][k] != -1) ) 
    { 
    adj[i][j] = adj[i][k] + adj[k][j]; 
    } 
    } 
    } 
    } 
    } 
    int main() { 
    int n=4 ; 
    // -1 in graph represents there is no direct edge between that two vertices 
    // weights in graph can be negative but not -1(as it represents INFINITY) 
    vector<vector<int>> graph = { 
        {0,9,-4,-1}, 
        {6,0,-1,2}, 
        {-1,5,0,-1}, 
        {-1,-1,1,0} 
    }; 
 
    floydWarshall(graph,n); 
     
    for(int i = 0; i<graph.size(); i++) { 
        for(int j = 0; j<graph.size(); j++) { 
            cout<<graph[i][j]<<" \t "; 
        } 
        cout<<endl; 
    } 
    return 0; 
    }

**INPUT:-**

    vector<vector<int>> graph = {
    {0, 9, -4, -1},
    {6, 0, -1, 2},
    {-1, 5, 0, -1},
    {-1, -1, 1, 0}
    };

**OUTPUT:-**

    0 	 8 	 -4  10 	
    6 	 0 	 2 	 2 	
    11 	 5 	 0 	 7 	
    7 	 1 	 1 	 0 	




  







  






