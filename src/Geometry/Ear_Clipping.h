#ifndef GRAPH_EAR_CLIPING
#define GRAPH_EAR_CLIPING

#include <vector>
#include <cmath>
#include <functional>

using namespace std;

namespace Graphics {

typedef vector<double> Point;   // x, y, z
typedef vector<double> triangle;   // x, y, z

inline bool isEar(Point& a, Point& b, Point& c, vector<Point>& polygon) {
    static std::function<double(Point&, Point&, Point&)> crossProduct =
        [](Point& a, Point& b, Point& c) {
        return
            (b[0] - a[0]) * (c[1] - b[1]) -
            (b[1] - a[1]) * (c[0] - b[0]);
    };

    // if not a convex vertice
    if (crossProduct(a, b, c) < 0)
        return false;

    for (Point p : polygon) {
        if ((p[0] == a[0] && p[1] == a[1]) ||
            (p[0] == b[0] && p[1] == b[1]) ||
            (p[0] == c[0] && p[1] == c[1]))
            continue;

        // if in trangle
        if (crossProduct(a, b, p) > 0 &&
            crossProduct(b, c, p) > 0 &&
            crossProduct(c, a, p) > 0
        )
            return false;
    }
    return true;
}

inline vector<triangle>& earClippingTriangulation(vector<Point>& polygon_, vector<triangle>& triangleSet) {
    vector<Point> polygon = polygon_;

    while (polygon.size() >= 3) {
        int n = polygon.size();
        
        for (int i = 0; i < n; i++) {
            int a = (i + n - 1) % n,
                c = (i + 1) % n;

            if (polygon.size() == 3 || 
                isEar(polygon[a], polygon[i], polygon[c], polygon)) 
            {
                triangleSet.push_back({
                    polygon[a][0], polygon[a][1], 0,
                    polygon[i][0], polygon[i][1], 0,
                    polygon[c][0], polygon[c][1], 0,
                });
                polygon.erase(polygon.begin() + i); 
                n--;
                i--; 
            }
        }
    }
    return triangleSet;
}

}
#endif // !GRAPH_EAR_CLIPING

