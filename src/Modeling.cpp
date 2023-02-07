#include "Modeling.h"

/*
 * 旋转体
 */
void Modeling::Rotator(Point& center, Point& axis, vector<Point>& f, int pointNum, double st, double ed) {
	double dAngle = (ed - st) / pointNum;
	Point p1, p2, p3, p4;
	Mat<double> rotateMat(3, 3), preRotateMat(3, 3);
	vector<double> direction(3), delta(3);

	for (int i = 0; i <= pointNum; i++) {
		// calculate rotate matrix
		double angle = st + dAngle * i;

		delta = { cos(angle), sin(angle), 0 };
		cross(direction, axis, delta),
		normalize(direction);

		double
			x = direction[0],
			y = direction[1],
			z = direction[2],
			sign = x * y > 0 ? -1 : 1,
			e = sqrt(1 - z * z),
			b = -x * z / e,
			d = -y * z / e,
			a = sign * abs(y / e),
			c = abs(x / e);

		rotateMat = {
			a, b, x,
			c, d, y,
			0, e, z
		};

		// generate
		if (i != 0) {
			for (int i = 1; i < f.size(); i++) {
				p1 = { f[i - 1][0], f[i - 1][1], 0};
				p2 = { f[i][0], f[i][1], 0};
				p3 = p1;
				p4 = p2;

				mul(p1, rotateMat, p1);
				mul(p2, rotateMat, p2);
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
		preRotateMat = rotateMat;
	}
}


/*
 * 平移体
 */
void Modeling::Translator(Point& st, Point& ed, vector<Point>& f) {
	// calculate rotate matrix
	Mat<double> rotateMat(3, 3);
	vector<double> direction(3);

	sub(direction, ed, st);
	normalize(direction);

	double
		x = direction[0],
		y = direction[1],
		z = direction[2],
		sign = x * y > 0 ? -1 : 1,
		e = sqrt(1 - z * z),
		b = -x * z / e,
		d = -y * z / e,
		a = sign * abs(y / e),
		c = abs(x / e);

	rotateMat = { 
		a, b, x,
		c, d, y,
		0, e, z
	};

	// generate
	Point stPoint(3), edPoint(3), preStPoint(3), preEdPoint(3);

	int fn = f.size();
	Point pt(3);

	for (int i = 0; i < fn; i++) {
		pt = { f[i][0], f[i][1], 0 };

		mul(pt, rotateMat, pt);
		add(stPoint, st, pt);
		add(edPoint, ed, pt);

		if (i != 0) {
			Triangle(stPoint,    preStPoint, edPoint);
			Triangle(preStPoint, preEdPoint, edPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}

void Modeling::Translator(vector<Point>& path, vector<Point>& f) {
	int n = path.size();
	Point p1 = path[0], p2;

	for (int i = 1; i < n; i++) {
		p2 = path[i];
		Translator(p1, p2, f);
		p1 = p2;
	}
}


/* --------------------------------
 *		平面图形
 * -------------------------------- */

/*
 * Triangle
 */
void Modeling::Triangle(Point& p1, Point& p2, Point& p3) {
	triangle tri(9);

	for (int i = 0; i < 3; i++) {
		tri[i] = p1[i];
		tri[i + 3] = p2[i];
		tri[i + 6] = p3[i];
	}

	Object.push_back(tri);
}

/* 
 * Rectangle
 */
void Modeling::Rectangle(Point& c, double X, double Y) {
	Point p1(3), p2(3), p3(3);
	Triangle(
		p1 = { c[0] + X / 2, c[1] + Y / 2, c[2] },
		p2 = { c[0] + X / 2, c[1] - Y / 2, c[2] },
		p3 = { c[0] - X / 2, c[1] + Y / 2, c[2] }
	);
	Triangle(
		p1 = { c[0] - X / 2, c[1] - Y / 2, c[2] },
		p2 = { c[0] + X / 2, c[1] - Y / 2, c[2] },
		p3 = { c[0] - X / 2, c[1] + Y / 2, c[2] }
	);
}

/*
 * Quadrangle
 */
void Modeling::Quadrangle(Point& p1, Point& p2, Point& p3, Point& p4) {
	Triangle(p1, p2, p3);
	Triangle(p1, p3, p4);
}

/* 
 * Polygon 
 */
// Convex Polygon 
void Modeling::ConvexPolygon(vector<Point>& p) {
	int n = p.size();

	for (int k = 1; k <= (n + 2) / 3; k++)
		for (int i = 0; i <= n - 2 * k; i += 2 * k)
			Triangle(p[i], p[i + k], p[(i + 2 * k) % n]);
}

void Modeling::Polygon(vector<Point>& p) {
	int n = p.size();
	Point p1(3), p2(3), p3(3);

	for (int k = 1; k <= (n + 2) / 3; k++)
		for (int i = 0; i <= n - 2 * k; i += 2 * k)
			Triangle(
				p1 = { p[i][0], p[i][1], 0},
				p2 = { p[i + k][0], p[i + k][1], 0},
				p3 = { p[(i + 2 * k) % n][0], p[(i + 2 * k) % n][1], 0}
			);
}

/* 
 * 画圆 
 */
void Modeling::Circle(Point& center, double r, int pointNum, double angleSt, double angleEd) {
	double dAngle = (angleEd - angleSt) / pointNum;
	Point ps(3), pe(3);

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
void Modeling::Surface(Mat<double>& z, double xs, double xe, double ys, double ye, Point* direct) {
	Point p(3), pl(3), pu(3), plu(3); 

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

			if (x == 0 || y == 0 ||
				z(x - 1, y) == HUGE_VAL ||  
				z(x, y - 1) == HUGE_VAL ||  
				z(x - 1, y - 1) == HUGE_VAL
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
void Modeling::Tetrahedron(Point& p1, Point& p2, Point& p3, Point& p4) {
	Triangle(p1, p2, p3);
	Triangle(p2, p3, p4);
	Triangle(p3, p4, p1);
	Triangle(p4, p1, p2);
}

/* 
 * 画矩体 
 */
void Modeling::Cuboid(Point& pMin, Point& pMax) {
	Point pMinTmp[3], pMaxTmp[3];
	for (int i = 0; i < 3; i++) {
		pMinTmp[i] = pMin; pMinTmp[i][i] = pMax[i];
		pMaxTmp[i] = pMax; pMaxTmp[i][i] = pMin[i];
	}

	for (int i = 0; i < 3; i++) {
		Quadrangle(pMin, pMinTmp[i], pMaxTmp[(i + 2) % 3], pMinTmp[(i + 1) % 3]);
		Quadrangle(pMax, pMaxTmp[i], pMinTmp[(i + 2) % 3], pMaxTmp[(i + 1) % 3]);
	}
}

void Modeling::Cuboid(Point& center, double X, double Y, double Z) {
	Point delta(3), pMax(3), pMin(3);
	delta = { X / 2, Y / 2, Z / 2 };

	add(pMax, center, delta);
	sub(pMin, center, delta);

	Cuboid(pMin, pMax);
}

/* 画圆台 * /
void Modeling::Frustum(Point& st, Point& ed, double Rst, double Red, int pointNum) {
	// 计算 Rotate Matrix
	Point direction, rotateAxis, rotateMat, zAxis(3), tmp;
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
	Point stPoint, edPoint, preStPoint, preEdPoint, deltaVector(3);
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
}*/

/* 画球 */
void Modeling::Sphere(Point& center, double r, int ThetaNum, int PhiNum, 
	double thetaSt, double thetaEd, 
	double phiSt, double phiEd
) {
	Point point(3), pointU(3), pointL(3), pointUL(3);
	double
		dTheta = (thetaEd - thetaSt) / ThetaNum,
		dPhi   = (phiEd - phiSt)     / PhiNum;

	for (int i = 1; i <= ThetaNum; i++) {
		double theta = thetaSt + i * dTheta;

		for (int j = 1; j <= PhiNum; j++) {
			double phi = phiSt + j * dPhi;

			point = {
				r * cos(phi) * cos(theta) + center[0],
				r * cos(phi) * sin(theta) + center[1],
				r * sin(phi) + center[2]
			};
			pointU = {
				r * cos(phi - dPhi) * cos(theta) + center[0],
				r * cos(phi - dPhi) * sin(theta) + center[1],
				r * sin(phi - dPhi) + center[2]
			};
			pointL = {
				r * cos(phi) * cos(theta - dTheta) + center[0],
				r * cos(phi) * sin(theta - dTheta) + center[1],
				r * sin(phi) + center[2]
			};
			pointUL = {
				r * cos(phi - dPhi) * cos(theta - dTheta) + center[0],
				r * cos(phi - dPhi) * sin(theta - dTheta) + center[1],
				r * sin(phi - dPhi) + center[2]
			};

			Triangle(point,  pointU, pointL);
			Triangle(pointL, pointU, pointUL);
		}
	}
}

/* --------------------------------
 *		Modifier
 * -------------------------------- */

void Modeling::Array(int count, double dx, double dy, double dz) {
	int n = Object.size();
	Point p1(3), p2(3), p3(3), delta(3);
	delta = { dx, dy, dz };

	for (int k = 1; k < count; k++) {
		for (int tri = 0; tri < n; tri++) {
			for (int dim = 0; dim < 3; dim++) {
				p1[dim] = Object[tri][dim];
				p2[dim] = Object[tri][3 + dim];
				p3[dim] = Object[tri][6 + dim];

				p1[dim] += k * delta[dim];
				p2[dim] += k * delta[dim];
				p3[dim] += k * delta[dim];
			}
			Triangle(p1, p2, p3);
		}
	}
}