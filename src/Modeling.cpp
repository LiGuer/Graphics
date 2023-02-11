#include "Modeling.h"

/*
 * 存储文件
 */
void Modeling::writeModel(const char* fileName) {
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

/*
 * 旋转体
 */
void Modeling::Rotator(Point& center, Point& axis, vector<Point>& f, int pointNum, double st, double ed, int isClosed) {
	double dAngle = (ed - st) / pointNum;
	Point p1(3), p2(3), p3(3), p4(3);
	Mat<double> rotateMat_0(3, 3), rotateMat(3, 3), preRotateMat(3, 3), firstRotateMat(3, 3);
	vector<double> direction(3), delta(3), direction_2(3);

	normalize(axis);

	// calculate the first rotate matrix
	double
		x = axis[0],
		y = axis[1],
		z = axis[2],
		sign = x * y > 0 ? -1 : 1,
		e = sqrt(1 - z * z),
		b = -x * z / e,
		d = -y * z / e,
		a = sign * abs(y / e),
		c = abs(x / e);

	if (e < 1e-4)
		E(rotateMat_0);
	else 
		rotateMat_0 = {
			a, b, x,
			c, d, y,
			0, e, z
		};

	for (int i = 0; i <= pointNum; i++) {
		double angle = st + dAngle * i;
		delta = { cos(angle), sin(angle), 0 };
		normalize(mul(direction, rotateMat_0, delta));
		normalize(cross(direction_2, axis, direction));

		// calculate rotate matrix
		rotateMat = {
			direction[0], axis[0], direction_2[0],
			direction[1], axis[1], direction_2[1],
			direction[2], axis[2], direction_2[2]
		};

		// generate
		if (i != 0) {
			for (int j = 1; j < f.size(); j++) {
				p3 = p1 = { f[j - 1][0], f[j - 1][1], 0 };
				p4 = p2 = { f[j][0],     f[j][1],     0 };

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
		if (i == 0)
			firstRotateMat = rotateMat;
	}

	// closed
	if (isClosed) {
		Point p1, p2, p3;
		vector<vector<double>> tris;

		Graphics::earClippingTriangulation(f, tris);

		int n = tris.size();
		
		for (int i = 0; i < n; i++) {
			mul(p1, rotateMat, p1 = { tris[i][0], tris[i][1], tris[i][2] });
			mul(p2, rotateMat, p2 = { tris[i][3], tris[i][4], tris[i][5] });
			mul(p3, rotateMat, p3 = { tris[i][6], tris[i][7], tris[i][8] });

			Object.push_back({
				p1[0] + center[0], p1[1] + center[1], p1[2] + center[2],
				p2[0] + center[0], p2[1] + center[1], p2[2] + center[2],
				p3[0] + center[0], p3[1] + center[1], p3[2] + center[2]
			});

			mul(p1, firstRotateMat, p1 = { tris[i][0], tris[i][1], tris[i][2] });
			mul(p2, firstRotateMat, p2 = { tris[i][3], tris[i][4], tris[i][5] });
			mul(p3, firstRotateMat, p3 = { tris[i][6], tris[i][7], tris[i][8] });

			Object.push_back({
				p1[0] + center[0], p1[1] + center[1], p1[2] + center[2],
				p2[0] + center[0], p2[1] + center[1], p2[2] + center[2],
				p3[0] + center[0], p3[1] + center[1], p3[2] + center[2]
			});
		}
	}
}


/*
 * 平移体
 */
void Modeling::Translator(Point& st, Point& ed, vector<Point>& f, int isClosed) {
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

	if (e < 1e-4)
		E(rotateMat);
	else
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

	// closed
	if (isClosed) {
		Point p1, p2, p3;
		vector<vector<double>> tris;

		Graphics::earClippingTriangulation(f, tris);

		int n = tris.size();

		for (int i = 0; i < n; i++) {
			mul(p1, rotateMat, p1 = { tris[i][0], tris[i][1], tris[i][2] });
			mul(p2, rotateMat, p2 = { tris[i][3], tris[i][4], tris[i][5] });
			mul(p3, rotateMat, p3 = { tris[i][6], tris[i][7], tris[i][8] });

			Object.push_back({
				p1[0] + st[0], p1[1] + st[1], p1[2] + st[2],
				p2[0] + st[0], p2[1] + st[1], p2[2] + st[2],
				p3[0] + st[0], p3[1] + st[1], p3[2] + st[2]
			});
			Object.push_back({
				p1[0] + ed[0], p1[1] + ed[1], p1[2] + ed[2],
				p2[0] + ed[0], p2[1] + ed[1], p2[2] + ed[2],
				p3[0] + ed[0], p3[1] + ed[1], p3[2] + ed[2]
			});
		}
	}
}

void Modeling::Translator(vector<Point>& path, vector<Point>& f, int isClosed) {
	int n = path.size();
	Point p1 = path[0], p2;

	for (int i = 1; i < n; i++) {
		p2 = path[i];

		if (isClosed && (i == 1 || i == n - 1))
			Translator(p1, p2, f, true);
		else
			Translator(p1, p2, f, false);

		p1 = p2;
	}
}

/* Rotator + Translator */
void Modeling::Rotator_Translator(
	Point& center, Point& axis, vector<Point>& f,
	vector<double>& direction_, double length,
	int pointNum, double st, double ed
) {
	double dAngle = (ed - st) / pointNum;
	Point p(3), p1(3), p2(3), p3(3), p4(3);
	Mat<double> rotateMat_0(3, 3), rotateMat(3, 3), preRotateMat(3, 3);
	vector<double> direction(3), delta(3), direction_2(3);

	normalize(axis);

	// calculate the first rotate matrix
	double
		x = axis[0],
		y = axis[1],
		z = axis[2],
		sign = x * y > 0 ? -1 : 1,
		e = sqrt(1 - z * z),
		b = -x * z / e,
		d = -y * z / e,
		a = sign * abs(y / e),
		c = abs(x / e);

	if (e < 1e-4)
		E(rotateMat_0);
	else
		rotateMat_0 = {
			a, b, x,
			c, d, y,
			0, e, z
	};

	for (int i = 0; i <= pointNum; i++) {
		double angle = st + dAngle * i;
		delta = { cos(angle), sin(angle), 0 };
		normalize(mul(direction, rotateMat_0, delta));
		normalize(cross(direction_2, axis, direction));

		// calculate rotate matrix
		rotateMat = {
			direction[0], axis[0], direction_2[0],
			direction[1], axis[1], direction_2[1],
			direction[2], axis[2], direction_2[2]
		};

		// generate
		if (i != 0) {
			for (int j = 1; j < f.size(); j++) {
				p3 = p1 = { f[j - 1][0], f[j - 1][1], 0 };
				p4 = p2 = { f[j][0],     f[j][1],     0 };

				mul(p1, rotateMat, p1);
				mul(p2, rotateMat, p2);
				mul(p3, preRotateMat, p3);
				mul(p4, preRotateMat, p4);

				add(p1, p1, center);
				add(p2, p2, center);
				add(p3, p3, center);
				add(p4, p4, center);

				add(p1, p1, mul(p, i / (double)pointNum * length, direction_));
				add(p2, p2, mul(p, i / (double)pointNum * length, direction_));
				add(p3, p3, mul(p, (i - 1) / (double)pointNum * length, direction_));
				add(p4, p4, mul(p, (i - 1) / (double)pointNum * length, direction_));

				Triangle(p1, p2, p3);
				Triangle(p4, p3, p2);
			}
		}
		preRotateMat = rotateMat;
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

void Modeling::Polygon(Point& c, vector<Point>& p) {
	vector<vector<double>> tris;

	Graphics::earClippingTriangulation(p, tris);
	addTriangleSet(c, tris);
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

void Modeling::Cuboid(Point& center, vector<double>& direction, double L, double W, double H) {
	vector<Point> f;
	vector<double> p(3), st(3), ed(3);
	f = {
		{-W / 2,-H / 2},
		{+W / 2,-H / 2},
		{+W / 2,+H / 2},
		{-W / 2,+H / 2},
		{-W / 2,-H / 2},
	};

	normalize(direction);

	add(ed, center, mul(p, +L / 2.0, direction));
	add(st, center, mul(p, -L / 2.0, direction));
	Translator(st, ed, f);
}

/* 画圆台 */
void Modeling::Frustum(Point& st, Point& ed, double Rst, double Red, int pointNum) {
	;
}

/* 画球 */
void Modeling::Sphere(Point& center, double r, int pointNum) {
	vector<double> st(3), ed(3);
	vector<int> N;
	vector<vector<double>> triangleSet;
	double more = r / pointNum * 3;

	Graphics::MarchingCubes([&](double x, double y, double z) {
		return r * r - (x * x + y * y + z * z); 
		},
		st = {-r - more,-r - more,-r - more },
		ed = { r + more, r + more, r + more },
		N  = { pointNum, pointNum, pointNum },
		triangleSet
	);

	addTriangleSet(center, triangleSet);
}

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
void Modeling::addTriangleSet(Point& center, vector<triangle>& tris) {
	int n = tris.size();

	for (int i = 0; i < n; i++) {
		Object.push_back({
			tris[i][0] + center[0], tris[i][1] + center[1], tris[i][2] + center[2],
			tris[i][3] + center[0], tris[i][4] + center[1], tris[i][5] + center[2],
			tris[i][6] + center[0], tris[i][7] + center[1], tris[i][8] + center[2],
			});
	}
}

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