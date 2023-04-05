#ifndef GEOMETRY_BEZIERCURVE_H
#define GEOMETRY_BEZIERCURVE_H

#include <vector>

using namespace std;

inline vector<double>& BezierCurve(vector<vector<double>>& points, double t, vector<double>& ans) {
	int n = points.size();
	vector<vector<double>> c = points;

	for (int i = 1; i < n; i++)
		for (int k = 0; k < n - i; k++)
			for (int d = 0; d < c[0].size(); d++)
				c[k][d] = c[k][d] * (1 - t) + c[k + 1][d] * t;

	return ans = c[0];
}

inline vector<vector<double>>& BezierCurve(vector<vector<double>>& points, int n, vector<vector<double>>& ans) {
	ans.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		BezierCurve(points, i / (double)n, ans[i]);
	}
	return ans;
}

#endif