#include "Intersect.h"

using namespace Matrix;
/*#############################################################################
* 
*						求交点
* 
##############################################################################*/

double Intersect::eps = 1e-5;

/*
 * 线段、线段交点
 */
bool OnSegments(Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	if (std::min(p1(1), p2(1)) <= p3(1) 
	&&  std::max(p1(1), p2(1)) >= p3(1)
	&&  std::min(p1(2), p2(2)) <= p3(2) 
	&&  std::max(p1(2), p2(2)) >= p3(2)
	)
		return true;
	return false;
}

bool Intersect::Segments(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4) {
	double 
		d1 = (p1(1) - p3(1)) * (p4(2) - p3(2)) - (p4(1) - p3(1)) * (p1(2) - p3(2)),
		d2 = (p2(1) - p3(1)) * (p4(2) - p3(2)) - (p4(1) - p3(1)) * (p2(2) - p3(2)),
		d3 = (p3(1) - p1(1)) * (p2(2) - p1(2)) - (p2(1) - p1(1)) * (p3(2) - p1(2)),
		d4 = (p4(1) - p1(1)) * (p2(2) - p1(2)) - (p2(1) - p1(1)) * (p4(2) - p1(2));

	if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0))
	&&  ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
		return true;

	else if (d1 == 0 && OnSegments(p3, p4, p1))
		return true;

	else if (d2 == 0 && OnSegments(p3, p4, p2))
		return true;

	else if (d3 == 0 && OnSegments(p1, p2, p3))
		return true;

	else if (d4 == 0 && OnSegments(p1, p2, p4))
		return true;

	return false;
}

/* 
 * 射线、平面交点
 */
double Intersect::RayPlane(Mat<>& raySt, Mat<>& ray, Mat<>& a, double b) {
	double t = dot(a, ray);
	if (t < eps) 
		return DBL_MAX;
	double d = (dot(a, raySt) - b) / t;
	return d > 0 ? d : DBL_MAX;
}

//3D
double Intersect::RayPlane(Mat<>& raySt, Mat<>& ray, double& A, double& B, double& C, double& D) {
	double t = A * ray[0] + B * ray[1] + C * ray[2];
	if (t < eps) 
		return DBL_MAX;
	double d = -(A * raySt[0] + B * raySt[1] + C * raySt[2] + D) / t;
	return d > 0 ? d : DBL_MAX;
}

/*
 * 射线、平面图案
 */
double Intersect::RayPlaneShape(Mat<>& raySt, Mat<>& ray, Mat<>& center, Mat<>& normal, Mat<>& one, bool(*f)(double, double)) {
	double
		d = RayPlane(raySt, ray, normal, dot(normal, center));
	if (d == DBL_MAX) 
		return DBL_MAX;

	static Mat<> delta, tmp;
	sub(delta, add(delta, raySt, mul(delta, d, ray)), center);
	cross_(tmp, delta, one);
	return f(dot(delta, one), (dot(tmp, normal) > 0 ? 1 : -1) * norm(tmp)) ? d : DBL_MAX;
}

/*
 * 射线、圆交点
 */
double Intersect::RayCircle(Mat<>& raySt, Mat<>& ray, Mat<>& center, double& R, Mat<>& normal) {
	double
		d = RayPlane(raySt, ray, normal, dot(normal, center));
	if (d == DBL_MAX) 
		return DBL_MAX;

	static Mat<> tmp;
	add(tmp, raySt, mul(tmp, d, ray));
	sub(tmp, tmp, center);
	return norm(tmp) <= R ? d : DBL_MAX;
}

/*
 * 射线、三角形交点
 */
double Intersect::RayTriangle(Mat<>& raySt, Mat<>& ray, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	static Mat<> edge[2], tmp, p, q;
	sub(edge[0], p2, p1);
	sub(edge[1], p3, p1);

	// p & a & tmp
	static double a, u, v;
	a = dot(cross_(p, ray, edge[1]), edge[0]);

	if (a > 0)
		sub(tmp, raySt, p1);
	else
		sub(tmp, p1, raySt), a = -a;

	if (a < 1e-4)
		return DBL_MAX;								//射线与三角面平行

	// u & q & v
	u = dot(p, tmp) / a;

	if (u < 0 || u > 1)
		return DBL_MAX;

	v = dot(cross_(q, tmp, edge[0]), ray) / a;
	return (v < 0 || u + v > 1) ? DBL_MAX : dot(q, edge[1]) / a;
}

/*
 * 射线、球面交点
 */
double Intersect::RaySphere(Mat<>& raySt, Mat<>& ray, Mat<>& center, double& R) {
	static Mat<> rayStCenter;
	sub(rayStCenter, raySt, center);

	double
		A = dot(ray, ray),
		B = 2 * dot(ray, rayStCenter),
		C = dot(rayStCenter, rayStCenter) - R * R,
		Delta = B * B - 4 * A * C;

	if (Delta < 0)
		return DBL_MAX;									//有无交点

	Delta = sqrt(Delta);
	return (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
}

/*
 * 射线、球面图案交点
 */
double Intersect::RaySphere(Mat<>& raySt, Mat<>& ray, Mat<>& center, double& R, bool(*f)(double, double)) {
	static Mat<> rayStCenter;
	sub(rayStCenter, raySt, center);

	double
		A = dot(ray, ray),
		B = 2 * dot(ray, rayStCenter),
		C = dot(rayStCenter, rayStCenter) - R * R,
		Delta = B * B - 4 * A * C;

	if (Delta < 0)
		return DBL_MAX;									//有无交点

	Delta = sqrt(Delta);

	if (f != NULL) {
		static double d; 
		static Mat<> delta;

		if ((d = (-B - Delta) / (2 * A)) > 1e-4) {
			sub(delta, add(delta, raySt, mul(delta, d, ray)), center);
			normalize(delta);
			if (f(
				acos(delta[2]), 
				atan(delta[1] / delta[0]) + (delta[1] >= 0 ? PI / 2 : PI / 2 * 3)
			))
				return d;
		}
		if ((d = (-B + Delta) / (2 * A)) > 1e-4) {
			sub(delta, add(delta, raySt, mul(delta, d, ray)), center);
			normalize(delta);
			if (f(
				acos(delta[2]), 
				atan(delta[1] / delta[0]) + (delta[1] >= 0 ? PI / 2 : PI / 2 * 3)
			))
				return d;
		}
		return DBL_MAX;
	}
	return (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
}

/*
 * 射线、二次曲面交点
 */
double Intersect::RayQuadric(Mat<>& raySt, Mat<>& ray, Mat<>& center, Mat<>& G) {
	static Mat<> rayStCenter, tmp;
	sub(rayStCenter, raySt, center);

	double
		A = dot(ray, mul(tmp, G, ray)),
		B = 2 * dot(ray, mul(tmp, G, rayStCenter)),
		C = dot(rayStCenter, mul(tmp, G, rayStCenter)) - 1,
		Delta = B * B - 4 * A * C;

	if (Delta < 0)
		return DBL_MAX;									//有无交点

	Delta = sqrt(Delta);
	return (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
}

/*
 * 射线、矩体交点
 */
double Intersect::RayCuboid(Mat<>& raySt, Mat<>& ray, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	return DBL_MAX;
}

double Intersect::RayCuboid(Mat<>& raySt, Mat<>& ray, Mat<>& pmin, Mat<>& pmax) {
	double t0 = -DBL_MAX, t1 = DBL_MAX;

	for (int dim = 0; dim < 3; dim++) {
		if (fabs(ray[dim]) < eps && (raySt[dim] < pmin[dim] || raySt[dim] > pmax[dim])) {
			return DBL_MAX;
		}
		double
			t0t = (pmin[dim] - raySt[dim]) / ray[dim],
			t1t = (pmax[dim] - raySt[dim]) / ray[dim];
		if (t0t > t1t)
			std::swap(t0t, t1t);

		t0 = std::max(t0, t0t);
		t1 = std::min(t1, t1t);

		if (t0 > t1 || t1 < 0)
			return DBL_MAX;
	}
	return t0 >= 0 ? t0 : t1;
}

/*
 * 射线、圆环交点
 */
double Intersect::RayTorus(Mat<>& raySt, Mat<>& ray, Mat<>& center, double R, double r) {
	double dx = ray[0], dy = ray[1], dz = ray[2],
		x0 = raySt[0], y0 = raySt[1], z0 = raySt[2],
		a = 4 * R * R,
		b = 2 * (x0 * dx + y0 * dy + z0 * dz),
		c = (R * R - r * r) + (x0 * x0 + y0 * y0 + z0 * z0),
		A = 1,
		B = 2 * b,
		C = b * b + 2 * c - a * (dx * dx + dy * dy),
		D = 2 * b * c - a * 2 * (x0 * dx + y0 * dy),
		E = c * c - a * (x0 * x0 + y0 * y0);

	complex<double> coeff[5], roots[4];
	coeff[4] = A;
	coeff[3] = B;
	coeff[2] = C;
	coeff[1] = D;
	coeff[0] = E;

	Function::solveQuarticEquation(coeff, roots);

	double minn = DBL_MAX;
	for (int i = 0; i < 4; i++) {
		if (abs(roots[i].imag()) < 10e-4 && roots[i].real() > 10e-4) {
			minn = min(minn, roots->real());
		}
	}
	return minn;
}
/*
double Intersect::RayPolynomialSurface(Mat<>& raySt, Mat<>& ray, Tensor<>& A) {
	int n = A.dimNum, m = A.dim[0] - 1;

	double dx = ray[0], dy = ray[1], dz = ray[2],
		x0 = raySt[0], y0 = raySt[1], z0 = raySt[2];

	complex<double>* coeff, roots;
	coeff = (complex<double>*) calloc((n + 1) * sizeof(complex<double>));
	roots = (complex<double>*) calloc( n      * sizeof(complex<double>));

	while (1) {

	}

	if(n == 4)
		solveQuartic(coeff, roots);

	double minn = DBL_MAX;
	for (int i = 0; i < 4; i++) {
		if (abs(roots[i].imag()) < 10e-4 && roots[i].real() > 10e-4) {
			minn = min(minn, roots->real());
		}
	}
	return minn;
}*/