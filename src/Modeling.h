#ifndef MODELING_H
#define MODELING_H

#include <vector>
#include "../../../Math/src/Matrix/Matrix.h"
#include "GraphicsIO.h"
#include "MarchingCubes.h"

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

/*
 * 存储文件
 */
inline void writeModel(const char* fileName) {
	char head[80] = { 0 };
	Mat<float> p[3], fv;
	Mat<double> t;

	p[0].alloc(3, Object.size());
	p[1].alloc(3, Object.size());
	p[2].alloc(3, Object.size());

	t.alloc(3, Object.size()).fill(1);
	normalize(t);

	fv.alloc(3, Object.size());
	for (int i = 0; i < t.size(); i++)
		fv(i) = t(i);

	Mat<short> attr(Object.size());

	for (int tri = 0; tri < Object.size(); tri++)
		for (int poi = 0; poi < 3; poi++)
			for (int dim = 0; dim < 3; dim++)
				p[poi](dim, tri) = Object[tri][poi * 3 + dim];

	GraphicsIO::stlWrite(fileName, head, fv, p[0], p[1], p[2], attr);
}

/* 图形 */
void Rotator	(Point& center, Point& axis, vector<Point>& f, int pointNum, double st = 0, double ed = 2 * PI);		// 旋转体 
void Translator	(Point& st, Point& ed, vector<Point>& f);						// 平移体 
void Translator	(vector<Point>& path, vector<Point>& f);

/* 2D Graph */
void Triangle   (Point& p1, Point& p2, Point& p3);
void Rectangle	(Point& c, double X, double Y);	
void Quadrangle	(Point& p1, Point& p2, Point& p3, Point& p4);
void ConvexPolygon(vector<Point>& p);
void Polygon    (vector<Point>& p);
void Circle		(Point& center, double r, int pointNum, double angleSt = 0, double angleEd = 2 * PI);
void Surface	(Mat<double>& z, double xs, double xe, double ys, double ye, Point* direct);

/* 3D Graph */
void Tetrahedron(Point& p1, Point& p2, Point& p3, Point& p4);
void Cuboid		(Point& pMin, Point& pMax);
void Cuboid		(Point& center, double X, double Y, double Z);
void Frustum	(Point& st, Point& ed, double Rst, double Red, int pointNum);							// 画圆台
void Sphere		(Point& center, double r, int ThetaNum, int PhiNum, 									// 画球
	double thetaSt = 0, double thetaEd = 2 * PI, double phiSt = -PI / 2, double phiEd = PI / 2);

/* Modifier */
void Array(int count, double dx, double dy, double dz);

};
#endif