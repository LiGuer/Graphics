#include "Intersect.h"

/*#############################################################################
* 
*						求交点
* 
##############################################################################*/

/* 射线、平面交点
 * d = - (A x_0 + B y_0 + C z_0 + D) / (A a + B b + C c)
 */
double RayPlane(Mat<>& RaySt, Mat<>& Ray, double& A, double& B, double& C, double& D) {
	double t = A * Ray[0] + B * Ray[1] + C * Ray[2];
	if (t == 0) return DBL_MAX;
	double d = -(A * RaySt[0] + B * RaySt[1] + C * RaySt[2] + D) / t;
	return d > 0 ? d : DBL_MAX;
}

double RayPlaneShape(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, Mat<>& normal, Mat<>& one, bool(*f)(double, double)) {
	double
		D = -dot(normal, Center),
		d = RayPlane(RaySt, Ray, normal[0], normal[1], normal[2], D);
	if (d == DBL_MAX) return DBL_MAX;
	static Mat<> delta, tmp;
	sub(delta, add(delta, RaySt, mul(delta, d, Ray)), Center);
	cross_(tmp, delta, one);
	return f(dot(delta, one), (dot(tmp, normal) > 0 ? 1 : -1) * norm(tmp)) ? d : DBL_MAX;
}

/* 射线、圆交点
 */
double RayCircle(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, double& R, Mat<>& normal) {
	double D = -(normal[0] * Center[0] + normal[1] * Center[1] + normal[2] * Center[2]),
		d = RayPlane(RaySt, Ray, normal[0], normal[1], normal[2], D);
	if (d == DBL_MAX) return DBL_MAX;

	static Mat<> tmp;
	add(tmp, RaySt, mul(tmp, d, Ray));
	sub(tmp, tmp, Center);
	return norm(tmp) <= R ? d : DBL_MAX;
}

/* 射线、三角形交点
 *
	O + t D = (1 - u - v)V0 + u V1 + v V2
	[ -D  V1-V0  V2-V0] [ t  u  v ]' = O - V0
	T = O - V0    E1 = V1 - V0    E2 = V2 - V0
	[ -D  E1  E2 ] [ t  u  v ]' = T
	t = | T  E1  E2| / |-D  E1  E2|
	u = |-D   T  E2| / |-D  E1  E2|
	v = |-D  E1  E2| / |-D  E1  E2|
 */
double RayTriangle(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	static Mat<> edge[2], tmp, p, q;
	sub(edge[0], p2, p1);
	sub(edge[1], p3, p1);

	// p & a & tmp
	static double a, u, v;
	a = dot(cross_(p, Ray, edge[1]), edge[0]);

	if (a > 0)
		sub(tmp, RaySt, p1);
	else
		sub(tmp, p1, RaySt), a = -a;

	if (a < 1e-4)
		return DBL_MAX;								//射线与三角面平行

	// u & q & v
	u = dot(p, tmp) / a;

	if (u < 0 || u > 1)
		return DBL_MAX;

	v = dot(cross_(q, tmp, edge[0]), Ray) / a;
	return (v < 0 || u + v > 1) ? DBL_MAX : dot(q, edge[1]) / a;
}

/* 射线、球面交点
 * K = ( -b ± sqrt(Δ) ) / 2 a
   Δ = b^2 - 4ac = 4(Al ΔX + Bl ΔY + Cl ΔZ)^2 - 4(Al^2 + Bl^2 + Cl^2)(ΔX^2 + ΔY^2 + ΔZ^2 - R^2)
   若Δ≥0 有交点.
 */
double RaySphere(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R, bool(*f)(double, double)) {
	static Mat<> RayStCenter;
	sub(RayStCenter, RaySt, center);

	double
		A = dot(Ray, Ray),
		B = 2 * dot(Ray, RayStCenter),
		Delta = B * B - 4 * A * (dot(RayStCenter, RayStCenter) - R * R);

	if (Delta < 0)
		return DBL_MAX;									//有无交点

	Delta = sqrt(Delta);

	if (f != NULL) {
		static double d; static Mat<> delta;
		if ((d = (-B - Delta) / (2 * A)) > 1e-4) {
			sub(delta, add(delta, RaySt, mul(delta, d, Ray)), center);
			normalize(delta);
			if (f(acos(delta[2]), atan(delta[1] / delta[0]) + (delta[1] >= 0 ? PI / 2 : PI / 2 * 3)))
				return d;
		}
		if ((d = (-B + Delta) / (2 * A)) > 1e-4) {
			sub(delta, add(delta, RaySt, mul(delta, d, Ray)), center);
			normalize(delta);
			if (f(acos(delta[2]), atan(delta[1] / delta[0]) + (delta[1] >= 0 ? PI / 2 : PI / 2 * 3)))
				return d;
		}
		return DBL_MAX;
	}
	return (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
}

/* 射线、矩体交点
 * d = - (A x_0 + B y_0 + C z_0 + D) / (A a + B b + C c)
 */
double RayCuboid(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	return DBL_MAX;
}

double RayCuboid(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax) {
	double t0 = -DBL_MAX, t1 = DBL_MAX;
	for (int dim = 0; dim < 3; dim++) {
		if (fabs(Ray[dim]) < EPS && (RaySt[dim] < pmin[dim] || RaySt[dim] > pmax[dim])) {
			return DBL_MAX;
		}
		double
			t0t = (pmin[dim] - RaySt[dim]) / Ray[dim],
			t1t = (pmax[dim] - RaySt[dim]) / Ray[dim];
		if (t0t > t1t)
			std::swap(t0t, t1t);

		t0 = std::max(t0, t0t);
		t1 = std::min(t1, t1t);

		if (t0 > t1 || t1 < 0)
			return DBL_MAX;
	}
	return t0 >= 0 ? t0 : t1;
}
