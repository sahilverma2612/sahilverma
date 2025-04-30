#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

// Structure to represent a point in 2D
struct Point {
    int x, y;
};

// Function to find the side of a point with respect to a line
int findSide(Point p1, Point p2, Point p) {
    int val = (p.y - p1.y) * (p2.x - p1.x) -
              (p2.y - p1.y) * (p.x - p1.x);
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}

// Function to calculate distance between point p and line p1p2
int lineDist(Point p1, Point p2, Point p) {
    return abs((p.y - p1.y) * (p2.x - p1.x) -
               (p2.y - p1.y) * (p.x - p1.x));
}

// Recursive function to find points on the convex hull
void quickHull(vector<Point>& points, Point p1, Point p2, int side, vector<Point>& hull) {
    int index = -1;
    int maxDist = 0;

    for (int i = 0; i < points.size(); i++) {
        int temp = lineDist(p1, p2, points[i]);
        if (findSide(p1, p2, points[i]) == side && temp > maxDist) {
            index = i;
            maxDist = temp;
        }
    }

    if (index == -1) {
        hull.push_back(p1);
        hull.push_back(p2);
        return;
    }

    // Recur for two parts
    quickHull(points, points[index], p1, -findSide(points[index], p1, p2), hull);
    quickHull(points, points[index], p2, -findSide(points[index], p2, p1), hull);
}

void printHull(vector<Point>& points) {
    if (points.size() < 3) {
        cout << "Convex hull not possible\n";
        return;
    }

    int min_x = 0, max_x = 0;
    for (int i = 1; i < points.size(); i++) {
        if (points[i].x < points[min_x].x) min_x = i;
        if (points[i].x > points[max_x].x) max_x = i;
    }

    vector<Point> hull;

    // Recursively find convex hull on both sides
    quickHull(points, points[min_x], points[max_x], 1, hull);
    quickHull(points, points[min_x], points[max_x], -1, hull);

    // Remove duplicates
    sort(hull.begin(), hull.end(), [](Point a, Point b) {
        return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    });
    hull.erase(unique(hull.begin(), hull.end(), [](Point a, Point b) {
        return a.x == b.x && a.y == b.y;
    }), hull.end());

    cout << "Points in Convex Hull:\n";
    for (Point p : hull)
        cout << "(" << p.x << ", " << p.y << ")\n";
}

int main() {
    vector<Point> points = {{0, 3}, {1, 1}, {2, 2}, {4, 4},
                            {0, 0}, {1, 2}, {3, 1}, {3, 3}};
    printHull(points);
    return 0;
}
