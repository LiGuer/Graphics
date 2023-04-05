#ifndef MODELING_H
#define MODELING_H

#include <vector>
#include "../../../Math/src/Matrix/Matrix.h"
#include "GraphicsIO.h"
#include "./Geometry/MarchingCubes.h"
#include "./Geometry/Ear_Clipping.h"

#define PI 3.141592653589

using namespace std;
using namespace Matrix;

class Modeling {
public:

	typedef vector<double> Point;   // x, y, z
	typedef vector<double> triangle;  // p1, p2, p3

	vector<triangle> Object;

	inline Modeling& operator= (Modeling& a){
		Object = a.Object;
		return *this;
	}

/* 存储文件 */
void writeModel(const char* fileName);

/* 图形 */
void Rotator	(Point& center, Point& axis, vector<Point>& f, 
				 int pointNum, double st = 0, double ed = 2 * PI, int isClosed = false);		// 旋转体 
void Translator	(Point& st, Point& ed, vector<Point>& f, int isClosed = true);			// 平移体 
void Translator	(vector<Point>& path, vector<Point>& f, int isClosed = true);
void Rotator_Translator
                (Point& center, Point& axis, vector<Point>& f,
	             vector<double>& direction, double length,
				 int pointNum, double st = 0, double ed = 2 * PI);
/* 2D Points */
static vector<Point>& Circle(vector<Point>& points, double R, int N, int clockwise = -1) {
	for (int i = 0; i <= N; i++) {
		double angle = clockwise * i / (double) N * 2 * PI;
		points.push_back({ R * sin(angle), R * cos(angle) });
	}
	return points;
}

/* 2D Graph */
void Triangle   (Point& p1, Point& p2, Point& p3);
void Rectangle	(Point& c, double X, double Y);	
void Quadrangle	(Point& p1, Point& p2, Point& p3, Point& p4);
void ConvexPolygon(vector<Point>& p);
void Polygon    (Point& c, vector<Point>& p);
void Circle		(Point& center, double r, int pointNum, double angleSt = 0, double angleEd = 2 * PI);
void Surface	(Mat<double>& z, double xs, double xe, double ys, double ye, Point* direct);

/* 3D Graph */
void Tetrahedron(Point& p1, Point& p2, Point& p3, Point& p4);
void Cuboid		(Point& pMin, Point& pMax);
void Cuboid		(Point& center, double X, double Y, double Z);
void Cuboid		(Point& center, vector<double>& direction, double L, double W, double H);
void Frustum	(Point& st, Point& ed, double Rst, double Red, int pointNum);							// 画圆台
void Sphere		(Point& center, double r, int pointNum);
void Sphere		(Point& center, double r, int ThetaNum, int PhiNum, 
				 double thetaSt = 0, double thetaEd = 2 * PI, 
			     double phiSt = -PI / 2, double phiEd = PI / 2);

/* Modifier */
void addTriangleSet(Point& center, vector<triangle>& tris);
void Array(int count, double dx, double dy, double dz);

};
#endif