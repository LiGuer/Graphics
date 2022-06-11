#ifndef MODELING_H
#define MODELING_H

#include <vector>
#include "../../../LiGu_Math/src/Math/Matrix/Matrix.h"
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

/*
 * 旋转体 
 */
void Rotator(Mat<>& center, Mat<>& axis, Mat<>& f, int pointNum, double st = 0, double ed = 2 * PI) {
	//Rotate f
	Mat<> p1(3), p2(3), p3(3), p4(3),
		RotateMat, preRotateMat, RotateMat0, RotateMatTmp,
		rotateAxis, fAxis(3), tmp;
	fAxis = { 0, 1, 0 };

	if (axis[0] != 0 || axis[2] != 0) {
		rotate(
			cross(rotateAxis, axis, fAxis),
			-acos(dot(axis, fAxis) / norm(axis)),
			E(RotateMatTmp.alloc(4, 4))
		);
		block(RotateMat0, RotateMatTmp, 1, 3, 1, 3);
	}
	else E(RotateMat0.alloc(3, 3));

	//main
	double dAngle = (ed - st) / pointNum;
	
	for (int i = 0; i <= pointNum; i++) {
		// 计算 Rotate Matrix
		rotate(
			axis, st + i * dAngle,
			E(RotateMatTmp.alloc(4, 4))
		);
		block(RotateMat, RotateMatTmp, 1, 3, 1, 3);
		mul(RotateMat, RotateMat, RotateMat0);

		// 画旋转体
		if (i != 0) {
			for (int i = 1; i < f.cols; i++) {
				p1 = { f(0, i - 1), f(1, i - 1), 0 };
				p2 = { f(0, i), f(1, i), 0 };
				p3 = p1;
				p4 = p2;

				mul(p1, RotateMat, p1);
				mul(p2, RotateMat, p2);
				mul(p3, preRotateMat, p3);
				mul(p4, preRotateMat, p4);

				add(p1, p1, center);
				add(p2, p2, center);
				add(p3, p3, center);
				add(p4, p4, center);

				Triangle(p1, p2, p3);
				Triangle(p4, p3, p2);
			}
		}
		preRotateMat = RotateMat;
	}
}

/* 
 * 平移体 
 */
void Translator(Mat<>& st, Mat<>& ed, Mat<>& f) {
	// 计算 Rotate Matrix
	Mat<> direction, rotateAxis, rotateMat, zAxis(3), tmp;
	zAxis = { 0, 0, 1 };
	sub(direction, ed, st);

	if (direction[0] != 0 || direction[1] != 0) {
		rotate(
			cross(rotateAxis, direction, zAxis),
			-acos(dot(direction, zAxis) / norm(direction)),
			E(rotateMat.alloc(4, 4))
		);
		block(rotateMat, rotateMat, 1, 3, 1, 3);
	}
	else E(rotateMat.alloc(3, 3));

	// 画圆台
	Mat<> stPoint, edPoint, preStPoint, preEdPoint;
	tmp.alloc(3);

	for (int i = 0; i < f.cols; i++) {
		mul(tmp, rotateMat, tmp = { f(0, i), f(1, i), 0 });
		add(stPoint, st, tmp);
		add(edPoint, ed, tmp);
		if (i != 0) {
			Triangle(stPoint,	preStPoint, edPoint);
			Triangle(preStPoint,preEdPoint, edPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}

void Translator(Mat<>& path, Mat<>& f) {
	Mat<> p1, p2;
	block(p1, path, 0, path.rows, 0, 0);
	p2 = p1;

	for (int i = 0; i < path.cols; i++) {
		block(p2, path, 0, path.rows, i, i);
		Translator(p1, p2, f);
		p1 = p2;
	}
}

/* --------------------------------
 *		平面图形
 * -------------------------------- */

/* 
 * 画矩形
 */
void Rectangle(Mat<>& c, double X, double Y) {
	Mat<> p1(3), p2(3), p3(3);
	Triangle(
		p1 = { c(0) + X / 2, c(1) + Y / 2, c(2) },
		p2 = { c(0) + X / 2, c(1) - Y / 2, c(2) },
		p3 = { c(0) - X / 2, c(1) + Y / 2, c(2) }
	);
	Triangle(
		p1 = { c(0) - X / 2, c(1) - Y / 2, c(2) },
		p2 = { c(0) + X / 2, c(1) - Y / 2, c(2) },
		p3 = { c(0) - X / 2, c(1) + Y / 2, c(2) }
	);
}

/* 
 * 画四边形 
 */
void Quadrangle(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4) {
	Triangle(p1, p2, p3);
	Triangle(p1, p3, p4);
}

/* 
 * 画多边形 
 */
void Polygon(Mat<> p[], int n) {
	for (int k = 1; k <= (n + 2) / 3; k++)
		for (int i = 0; i <= n - 2 * k; i += 2 * k)
			Triangle(p[i], p[i + k], p[(i + 2 * k) % n]);
}

/* 
 * 画圆 
 */
void Circle(Mat<>& center, double r, int pointNum, double angleSt = 0, double angleEd = 2 * PI) {
	double dAngle = (angleEd - angleSt) / pointNum;
	Mat<> ps(3), pe(3);

	for (int i = 0; i < pointNum; i++) {
		double theta = i * dAngle;
		ps = {
			r * cos(theta + angleSt),
			r * sin(theta + angleSt),
			0
		};
		pe = {
			r * cos(theta + angleSt + dAngle),
			r * sin(theta + angleSt + dAngle),
			0
		};
		Triangle(
			add(ps, ps, center),
			add(pe, pe, center),
			center
		);
	}
}

/* 
 * 画曲面 
 */
void Surface(Mat<>& z, double xs, double xe, double ys, double ye, Mat<>* direct) {
	Mat<> p(3), pl(3), pu(3), plu(3), FaceVec, tmp, light(3); 
	light = 1 / sqrt(3);

	double 
		dx = (xe - xs) / z.rows,
		dy = (ye - ys) / z.cols;

	for (int y = 0; y < z.cols; y++) {
		for (int x = 0; x < z.rows; x++) {
			if (z(x, y) == HUGE_VAL) 
				continue;

			p = { 
				xs + x * dx, 
				ys + y * dy, 
				z(x, y)
			};

			if (x == 0 || y == 0) 
				continue;

			if (z(x - 1, y) == HUGE_VAL
			||  z(x, y - 1) == HUGE_VAL
			||  z(x - 1, y - 1) == HUGE_VAL
			) 
				continue;

			pl = { 
				xs + (x - 1) * dx, 
				ys + y * dy,
				z(x - 1, y) 
			};
			pu = {
				xs + x * dx,
				ys + (y - 1) * dy, 
				z(x, y - 1) 
			};
			plu = { 
				xs + (x - 1) * dx,	
				ys + (y - 1) * dy,
				z(x - 1, y - 1) 
			};

			Triangle(p, pl, pu);
			Triangle(plu, pu, pl);
		}
	}
}

/* --------------------------------
 *		三维图形
 * -------------------------------- */

/* 
 * 画四面体 
 */
void Tetrahedron(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4) {
	Triangle(p1, p2, p3);
	Triangle(p2, p3, p4);
	Triangle(p3, p4, p1);
	Triangle(p4, p1, p2);
}

/* 
 * 画矩体 
 */
void Cuboid(Mat<>& pMin, Mat<>& pMax) {
	Mat<> pMinTmp[3], pMaxTmp[3];
	for (int i = 0; i < 3; i++) {
		pMinTmp[i] = pMin; pMinTmp[i][i] = pMax[i];
		pMaxTmp[i] = pMax; pMaxTmp[i][i] = pMin[i];
	}

	for (int i = 0; i < 3; i++) {
		Quadrangle(pMin, pMinTmp[i], pMaxTmp[(i + 2) % 3], pMinTmp[(i + 1) % 3]);
		Quadrangle(pMax, pMaxTmp[i], pMinTmp[(i + 2) % 3], pMaxTmp[(i + 1) % 3]);
	}
}

void Cuboid(Mat<>& center, double X, double Y, double Z) {
	Mat<> delta(3), pMax(3), pMin(3);
	delta = { X, Y, Z };

	add(pMax, center, delta);
	sub(pMin, center, delta);

	Cuboid(pMin, pMax);
}

/* 画圆台 */
void Frustum(Mat<>& st, Mat<>& ed, double Rst, double Red, int pointNum) {
	// 计算 Rotate Matrix
	Mat<> direction, rotateAxis, rotateMat, zAxis(3), tmp;
	zAxis = { 0, 0, 1 };
	sub(direction, ed, st);

	if (direction[0] != 0 || direction[1] != 0) {
		rotate(
			cross(rotateAxis, direction, zAxis),
			-acos(dot(direction, zAxis) / norm(direction)),
			E(rotateMat.alloc(4, 4))
		);
		block(rotateMat, rotateMat, 1, 3, 1, 3);
	}
	else E(rotateMat.alloc(3, 3));

	// 画圆台
	Mat<> stPoint, edPoint, preStPoint, preEdPoint, deltaVector(3);
	double dAngle = 2.0 * PI / pointNum;

	for (int i = 0; i <= pointNum; i++) {
		deltaVector = {
			cos(i * dAngle),
			sin(i * dAngle),
			0
		};

		mul(deltaVector, rotateMat, deltaVector);
		add(stPoint, st, mul(stPoint, Rst, deltaVector));
		add(edPoint, ed, mul(edPoint, Red, deltaVector));

		if (i != 0) {
			Triangle(stPoint, preStPoint, edPoint);
			Triangle(preStPoint, preEdPoint, edPoint);
			Triangle(st, stPoint, preStPoint);
			Triangle(ed, edPoint, preEdPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}

/* 画球 */
void Sphere(Mat<>& center, double r, double dAngle,
	double thetaSt = 0, double thetaEd = 2 * PI, 
	double phiSt = -PI / 2, double phiEd = PI / 2
) {
	Mat<> point(3), pointU(3), pointL(3), pointUL(3);
	int ThetaNum = (thetaEd - thetaSt) / dAngle,
		PhiNum = (phiEd - phiSt) / dAngle;

	for (int i = 1; i <= ThetaNum; i++) {
		double theta = thetaSt + i * dAngle;

		for (int j = 1; j <= PhiNum; j++) {
			double phi = phiSt + j * dAngle;

			point = {
				r * cos(phi) * cos(theta) + center[0],
				r * cos(phi) * sin(theta) + center[1],
				r * sin(phi) + center[2]
			};
			pointU = {
				r * cos(phi - dAngle) * cos(theta) + center[0],
				r * cos(phi - dAngle) * sin(theta) + center[1],
				r * sin(phi - dAngle) + center[2]
			};
			pointL = {
				r * cos(phi) * cos(theta - dAngle) + center[0],
				r * cos(phi) * sin(theta - dAngle) + center[1],
				r * sin(phi) + center[2]
			};
			pointUL = {
				r * cos(phi - dAngle) * cos(theta - dAngle) + center[0],
				r * cos(phi - dAngle) * sin(theta - dAngle) + center[1],
				r * sin(phi - dAngle) + center[2]
			};

			Triangle(point,  pointU, pointL);
			Triangle(pointL, pointU, pointUL);
		}
	}
}

/* --------------------------------
 *		Modifier
 * -------------------------------- */

void Array(int count, double dx, double dy, double dz) {
	int n = size();
	Mat<> p1(3), p2(3), p3(3), delta(3);
	delta = { dx, dy, dz };

	for (int k = 1; k < count; k++) {
		for (int i = 0; i < n; i++) {
			for (int dim = 0; dim < 3; dim++) {
				p1(dim) = Object[i * 9 + dim];
				p2(dim) = Object[i * 9 + 3 + dim];
				p3(dim) = Object[i * 9 + 6 + dim];

				p1(dim) += k * delta(dim);
				p2(dim) += k * delta(dim);
				p3(dim) += k * delta(dim);
			}
			Triangle(p1, p2, p3);
		}
	}
}


};
#endif