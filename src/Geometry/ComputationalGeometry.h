#ifndef COMPUTATIONAL_GEOMETRY_H
#define COMPUTATIONAL_GEOMETRY_H

#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <vector>
#include"../Matrix/Matrix.h"

using namespace Matrix;

#define PI 3.141592653589
#define EPS 10e-4

namespace Geometry {

/****************************************************************
				是否在三角内
****************************************************************/
inline bool inTriangle(Mat<>& p0, Mat<>& TriP1, Mat<>& TriP2, Mat<>& TriP3) {
	Mat<> tmp, edge[2];
	sub(edge[0],TriP2, TriP1);
	sub(edge[1],TriP3, TriP1);
	sub(tmp,       p0, TriP1);

	double 
		Dot00 = dot(edge[0], edge[0]),
		Dot01 = dot(edge[0], edge[1]),
		Dot11 = dot(edge[1], edge[1]),
		Dot02 = dot(edge[0], tmp),
		Dot12 = dot(edge[1], tmp),
		t = Dot00 * Dot11 - Dot01 * Dot01,
		u =(Dot11 * Dot02 - Dot01 * Dot12) / t,
	    v =(Dot00 * Dot12 - Dot01 * Dot02) / t;
	return (u < 0 || u > 1 || v < 0 || v > 1 || u + v > 1) ? false : true;
}

/****************************************************************
				判断四点共圆
*	[输出]: 圆外-1，圆上0，圆内1
*	三点确定圆方程: 即 解行列式:
		| x1²+y1²  x1  y1  1 | ?= 0
		| x2²+y2²  x2  y2  1 |
		| x3²+y3²  x3  y3  1 |
		| x4²+y4²  x4  y4  1 |
*	[几何解释]: 通过把平面点提升到三维的抛物面中，由于抛物面被平面所截的截面为圆形，四点共面即使共圆，也可以用四面体的体积判断是否共圆。
****************************************************************/
inline bool onCircle(Mat<> Points[]) {
	Mat<> mat(4, 4);
	for (int i = 0; i < 4; i++) {
		mat(i, 0) = dot(Points[i], Points[i]);
		mat(i, 1) = Points[i][0];
		mat(i, 2) = Points[i][1];
		mat(i, 4) = 1;
	}
	return det(mat) < EPS ? true : false;
}

/****************************************************************
				平面三点确定圆方程
*	[公式]: 圆方程: (x - cx)² + (y - cy)² = R²
*	[算法]: 三点确定圆方程: 即 解行列式:
			| x²+y²    x   y   1 |  =  0
			| x1²+y1²  x1  y1  1 |
			| x2²+y2²  x2  y2  1 |
			| x3²+y3²  x3  y3  1 |
		即.目标三点和圆上(x,y)应该满足方程组:
			(x²+y²)·a + x·b + y·c + 1·d = 0
*	[推导]:
			M11(x²+y²) - M12 x + M13 y - M14 = (x²+y²)·a + x·b + y·c + 1·d = 0
			=> (x² + b/a x) + (y² + c/a y) = - d/a
			=> (x + b/2a)² + (y + c/2a)² = -d/a + b²/4a² + c²/4a²
							CircumCircle 三角形外接圆
*	外接圆圆心: 即. 三点确定圆方程问题， 也是任意两边的垂直平分线的交点.直接用 ThreePointsToCircle()方法
****************************************************************/
inline Mat<>& ThreePoints2Circle(Mat<> Points[], Mat<>& center, double& R) {
	Mat<> mat(4, 4);
	for (int i = 0; i < 3; i++) {
		mat(i + 1, 0) = dot(Points[i], Points[i]);
		mat(i + 1, 1) = Points[i][0];
		mat(i + 1, 2) = Points[i][1];
		mat(i + 1, 3) = 1;
	}
	double 
		a =  comi(mat, 0, 0),
		b = -comi(mat, 0, 1),
		c =  comi(mat, 0, 2),
		d = -comi(mat, 0, 3);
	R = sqrt(-d / a + b * b / (4 * a * a) + c * c / (4 * a * a));

	center.zero(2);
	return center = { -b / (2 * a), -c / (2 * a) };
}

/****************************************************************
				球面均匀点分布
*	[Referance]:
		[1] Thanks and copyright for https://github.com/SebLague/Boids
****************************************************************/
inline Mat<>* getSphereFibonacciPoint(int n) {
	Mat<>* point = (Mat<>*)malloc(n * sizeof(Mat<>));
	memset(point, 0, n * sizeof(Mat<>));
	double goldenRatio = (1 + sqrt(5)) / 2, angleIncrement = PI * 2 * goldenRatio;	// 黄金分割点
	for (int i = 0; i < n; i++) {
		double t = (double)i / n, inclination = acos(1 - 2 * t), azimuth = angleIncrement * i;
		point[i].zero(3);
		point[i][0] = sin(inclination) * cos(azimuth);
		point[i][1] = sin(inclination) * sin(azimuth);
		point[i][2] = cos(inclination);
	} return point;
}

}
#endif