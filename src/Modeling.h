#ifndef MODELING_H
#define MODELING_H

#include <vector>
#include "../../../../../Math/src/Math/Matrix/Matrix.h"
#include "../../../../../Math/src/Math/Geometry/ComputationalGeometry.h"
#include "../GraphicsIO.h"

#define PI 3.141592653589

using namespace Matrix;

class Modeling {
public:

	std::vector<double> Object;

	inline int size() {
		return Object.size() / 9;
	}

	inline Modeling& operator= (Modeling& a){
		Object = a.Object;
		return *this;
	}

/*
 * 三角形
 */
void Triangle(Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	for(int i = 0; i < 3; i++)
		Object.push_back(p1(i));

	for(int i = 0; i < 3; i++)
		Object.push_back(p2(i));

	for(int i = 0; i < 3; i++)
		Object.push_back(p3(i));
}

/*
 * 存储文件
 */
void writeModel(const char* fileName) {
	char head[80] = { 0 };
	Mat<float> p[3], fv;
	Mat<> t;

	p[0].alloc(3, size());
	p[1].alloc(3, size());
	p[2].alloc(3, size());

	t.alloc(3, size()).fill(1);
	normalize(t);

	fv.alloc(3, size());
	for (int i = 0; i < t.size(); i++)
		fv(i) = t(i);

	Mat<short> attr(size());

	for (int i = 0; i < Object.size(); i++)
		p[(i % 9) / 3](i % 3, i / 9) = Object[i];

	GraphicsIO::stlWrite(fileName, head, fv, p[0], p[1], p[2], attr);
}

/* 图形 */
void Rotator	(Mat<>& center, Mat<>& axis, Mat<>& f, int pointNum, double st = 0, double ed = 2 * PI);		// 旋转体 
void Translator	(Mat<>& st, Mat<>& ed, Mat<>& f);						// 平移体 
void Translator	(Mat<>& path, Mat<>& f);

/* 平面图形 */
void Rectangle	(Mat<>& c, double X, double Y);							// 画矩形
void Quadrangle	(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4);			// 画四边形 
void ConvexPolygon(Mat<>* p, int n);									// 画凸多边形
void Polygon	(Mat<>* p, int n);										// 画多边形
void Polygon(Mat<>& p);
void Circle		(Mat<>& center, double r, int pointNum, double angleSt = 0, double angleEd = 2 * PI);	// 画圆 
void Surface	(Mat<>& z, double xs, double xe, double ys, double ye, Mat<>* direct);					// 画曲面

/* 三维图形 */
void Tetrahedron(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4);											// 画四面体
void Cuboid		(Mat<>& pMin, Mat<>& pMax);																// 画矩体 
void Cuboid		(Mat<>& center, double X, double Y, double Z);
void Frustum	(Mat<>& st, Mat<>& ed, double Rst, double Red, int pointNum);							// 画圆台
void Sphere		(Mat<>& center, double r, int ThetaNum, int PhiNum, 									// 画球
	double thetaSt = 0, double thetaEd = 2 * PI, double phiSt = -PI / 2, double phiEd = PI / 2);

/* Modifier */
void Array(int count, double dx, double dy, double dz);

};
#endif