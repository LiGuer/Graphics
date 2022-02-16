/*
Copyright 2020,2021 LiGuer. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
[Reference]:
[1] Thanks for Kevin Beason at http://www.kevinbeason.com/smallpt/
==============================================================================*/
#include "RayTracing.h"

/*#############################################################################
*						几何光学  
##############################################################################*/

/* 反射
 * 入射角 == 反射角 $\.L_o = \.L_i - 2 (\.N · \.L_i) \.N$
 */
Mat<>& reflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
	return RayO.add(RayO.mul(-2 * faceVec.dot(RayI), faceVec), RayI).normalize();
}

/* 折射
 * $n_i·\sin \theta_i = n_o·\sin \theta_o$
   $\.L_o = \.L_i + \.N· (\cos\theta_i - \frac{\sqrt{1 - n^2·\sin^2\theta_i}}{n})$
 */
Mat<>& refract(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO) {
	double k = rateI / rateO,
		CosI = faceVec.dot(RayI),
		CosO = 1 - pow(k, 2) * (1 - pow(CosI, 2));
	return CosO < 0 ? reflect(RayI, faceVec, RayO) :				//全反射
		RayO.add(RayI, RayO.mul(-CosI - (CosI > 0 ? -1 : 1)* sqrt(CosO) / k, faceVec)).normalize();
}

/* 漫反射
 * MonteCarlo: 随机持续采样, 在面矢半球内，面积均匀的随机取一射线，作为反射光线.
 */
Mat<>& diffuseReflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
	double r1 = 2 * PI * RAND_DBL, r2 = RAND_DBL;
	static Mat<> t(3), u, v;
	faceVec *= faceVec.dot(RayI) > 0 ? -1 : 1;
	t[0] = fabs(faceVec[0]) > 0.1 ? 0 : 1;
	t[1] = t[0] == 0 ? 1 : 0;
	u.mul(cos(r1) * sqrt(r2), u.cross_(t, faceVec).normalize());
	v.mul(sin(r1) * sqrt(r2), v.cross_(faceVec, u).normalize());
	return RayO.add(RayO.mul(sqrt(1 - r2), faceVec), u += v).normalize();
}

/* 雾 (均匀同质)
 * $L_o(x) = L_i(x) · t(x) + A · (1 - t(x))$
   $t(x) = e^{-\beta d(x)}$
 */
double Haze(double I, double A, double dis, double beta) {
	double t = exp(-beta * dis);
	return I * t + A * (1 - t);
}
Mat<>& Haze(Mat<>& I, Mat<>& O, Mat<>& A, double dis, double beta) {
	static Mat<> Tmp(3);
	double t = exp(-beta * dis);
	O.mul(t, I);
	Tmp.mul(1 - t, A);
	O += Tmp;
	return O;
}

/*#############################################################################
						求交点
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
		D = -normal.dot(Center),
		d = RayPlane(RaySt, Ray, normal[0], normal[1], normal[2], D);
	if (d == DBL_MAX) return DBL_MAX;
	static Mat<> delta, tmp; 
	delta.sub(delta.add(RaySt, delta.mul(d, Ray)), Center);
	tmp.cross_(delta, one);
	return f(delta.dot(one), (tmp.dot(normal) > 0 ? 1 : -1) * tmp.norm()) ? d : DBL_MAX;
}

/* 射线、圆交点
 */
double RayCircle(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, double& R, Mat<>& normal) {
	double D = -(normal[0] * Center[0] + normal[1] * Center[1] + normal[2] * Center[2]),
		   d = RayPlane(RaySt, Ray, normal[0], normal[1], normal[2], D);
	if (d == DBL_MAX) return DBL_MAX;
	static Mat<> tmp; tmp.add(RaySt, tmp.mul(d, Ray)); 
	return (tmp -= Center).norm() <= R ? d : DBL_MAX;
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
	edge[0].sub(p2, p1);
	edge[1].sub(p3, p1);
	// p & a & tmp
	static double a, u, v;
	a = p.cross_(Ray, edge[1]).dot(edge[0]);
	if (a > 0) tmp.sub(RaySt, p1);
	else       tmp.sub(p1, RaySt), a = -a;
	if (a < 1e-4)		return DBL_MAX;								//射线与三角面平行
	// u & q & v
	u = p.dot(tmp) / a;
	if (u < 0 || u > 1)	return DBL_MAX;
	v = q.cross_(tmp, edge[0]).dot(Ray) / a;
	return (v < 0 || u + v > 1) ? DBL_MAX : q.dot(edge[1]) / a;
}

/* 射线、球面交点
 * K = ( -b ± sqrt(Δ) ) / 2 a
   Δ = b^2 - 4ac = 4(Al ΔX + Bl ΔY + Cl ΔZ)^2 - 4(Al^2 + Bl^2 + Cl^2)(ΔX^2 + ΔY^2 + ΔZ^2 - R^2)
   若Δ≥0 有交点.
 */
double RaySphere(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R, bool(*f)(double, double)) {
	static Mat<> RayStCenter; RayStCenter.sub(RaySt, center);
	double 
		A = Ray.dot(Ray),
		B = 2 * Ray.dot(RayStCenter),
		Delta = B * B - 4 * A * (RayStCenter.dot(RayStCenter) - R * R);
	if (Delta < 0) return DBL_MAX;									//有无交点
	Delta = sqrt(Delta);
	if (f != NULL) {
		static double d; static Mat<> delta;
		if ((d = (-B - Delta) / (2 * A)) > 1e-4) {
			delta.sub(delta.add(RaySt, delta.mul(d, Ray)), center).normalize();
			if (f(acos(delta[2]), atan(delta[1] / delta[0]) + (delta[1] >= 0 ? PI / 2 : PI / 2 * 3))) return d;
		}
		if ((d = (-B + Delta) / (2 * A)) > 1e-4) {
			delta.sub(delta.add(RaySt, delta.mul(d, Ray)), center).normalize();
			if (f(acos(delta[2]), atan(delta[1] / delta[0]) + (delta[1] >= 0 ? PI / 2 : PI / 2 * 3))) return d;
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
		if (t0t > t1t) std::swap(t0t, t1t);
		t0 = std::max(t0, t0t);
		t1 = std::min(t1, t1t);
		if (t0 > t1 || t1 < 0) return DBL_MAX;
	}
	return t0 >= 0 ? t0 : t1;
}

/*#############################################################################
						对象/对象树
##############################################################################*/
/*--------------------------------[ 建树 ]--------------------------------*/
void ObjectTree::build(std::vector<Object>& obSet) {
	ObNodeList = (ObjectNode*)malloc(obSet.size() * sizeof(ObjectNode));
	Mat<> p(3);
	for (int i = 0; i < obSet.size(); i++) {
		Object* bound = new Object;
		bound->type = CUBOID;
		bound->v = (void**)calloc(2, sizeof(void*));
		bound->v[0] = new Mat<>(3);
		bound->v[1] = new Mat<>(3);

		Object* ob = &obSet[i];
		switch (ob->type) {
		case CIRCLE:
			p = *(double*)ob->v[2];
			(*(Mat<>*)bound->v[0]).sub(*(Mat<>*)ob->v[0], p);
			(*(Mat<>*)bound->v[1]).add(*(Mat<>*)ob->v[0], p);
			break;
		case TRIANGLE:
			for (int j = 0; j < 3; j++) {
				(*(Mat<>*)bound->v[0])[j] = std::min((*(Mat<>*)ob->v[0])[j], std::min((*(Mat<>*)ob->v[1])[j], (*(Mat<>*)ob->v[2])[j]));
				(*(Mat<>*)bound->v[1])[j] = std::max((*(Mat<>*)ob->v[0])[j], std::max((*(Mat<>*)ob->v[1])[j], (*(Mat<>*)ob->v[2])[j]));
			}
			break;
		case SPHERE:
			p = *(double*)ob->v[1];
			(*(Mat<>*)bound->v[0]).sub(*(Mat<>*)ob->v[0], p);
			(*(Mat<>*)bound->v[1]).add(*(Mat<>*)ob->v[0], p);
			break;
		case CUBOID: delete bound; bound = ob; break;
		}
		ObNodeList[i].ob = &obSet[i];
		ObNodeList[i].bound = bound;
	}
	std::sort(ObNodeList, ObNodeList + obSet.size(), [](ObjectNode& a, ObjectNode& b) {
		if (a.ob->type == PLANE || a.ob->type == PLANESHAPE) return true; return false;
	});
	while (ObNodeList[planeNum].ob->type == PLANE || ObNodeList[planeNum].ob->type == PLANESHAPE) planeNum++; 
	build(ObNodeList, planeNum, obSet.size() - 1, root); 
}

void ObjectTree::build(ObjectNode* obSet, int l, int r, ObjectNode*& node) {
	if (l == r) { node = &obSet[l]; return; }
	node = new ObjectNode;
	Object* bound = new Object;
	bound->type = CUBOID;
	bound->v = (void**)calloc(2, sizeof(void*));
	bound->v[0] = new Mat<>; *(Mat<>*)bound->v[0] = *(Mat<>*)obSet[l].bound->v[0];
	bound->v[1] = new Mat<>; *(Mat<>*)bound->v[1] = *(Mat<>*)obSet[l].bound->v[1];

	double delta[3];
	for (int i = l + 1; i <= r; i++) {
		for (int j = 0; j < 3; j++) {
			(*(Mat<>*)bound->v[0])[j] = std::min((*(Mat<>*)bound->v[0])[j], (*(Mat<>*)obSet[i].bound->v[0])[j]);
			(*(Mat<>*)bound->v[1])[j] = std::max((*(Mat<>*)bound->v[1])[j], (*(Mat<>*)obSet[i].bound->v[1])[j]);
			delta[j] = std::max(delta[j], (*(Mat<>*)obSet[i].bound->v[1])[j] - (*(Mat<>*)obSet[i].bound->v[0])[j]);
		}
	}
	node->bound = bound;
	int dim;
	dim = delta[0] / ((*(Mat<>*)bound->v[1])[0] - (*(Mat<>*)bound->v[0])[0]) < delta[1] / ((*(Mat<>*)bound->v[1])[1] - (*(Mat<>*)bound->v[0])[1]) ? 0 : 1;
	dim = delta[1] / ((*(Mat<>*)bound->v[1])[1] - (*(Mat<>*)bound->v[0])[1]) < delta[2] / ((*(Mat<>*)bound->v[1])[2] - (*(Mat<>*)bound->v[0])[2]) ? 1 : 2;
	std::sort(obSet + l, obSet + r + 1, [&dim](ObjectNode& a, ObjectNode& b) {
		if((*(Mat<>*)a.bound->v[0])[dim] != (*(Mat<>*)b.bound->v[0])[dim]) 
			return (*(Mat<>*)a.bound->v[0])[dim] < (*(Mat<>*)b.bound->v[0])[dim];
			return (*(Mat<>*)a.bound->v[1])[dim] < (*(Mat<>*)b.bound->v[1])[dim];
		});
	build(obSet, l, (l + r) / 2,     node->kid[0]);
	build(obSet, (l + r) / 2 + 1, r, node->kid[1]);
}

/*--------------------------------[ 求交 ]--------------------------------*/
double ObjectTree::seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object*& ob) {
	double disMin = seekIntersection(RaySt, Ray, root, ob), dis_t;
	for (int i = 0; i < planeNum; i++) {
		dis_t = seekIntersection(RaySt, Ray, *ObNodeList[i].ob);
		if (dis_t > EPS && disMin > dis_t) { ob = ObNodeList[i].ob; disMin = dis_t; }
	}
	return disMin;
}

double ObjectTree::seekIntersection(Mat<>& RaySt, Mat<>& Ray, ObjectNode* node, Object*& ob) {
	if (node->ob != NULL) { ob = node->ob; return seekIntersection(RaySt, Ray, *node->ob); }
	if (seekIntersection(RaySt, Ray, *node->bound) == DBL_MAX) return DBL_MAX;
	Object* ob_1, * ob_2;
	double dis_1 = seekIntersection(RaySt, Ray, node->kid[0], ob_1); dis_1 = dis_1 > EPS ? dis_1 : DBL_MAX;
	double dis_2 = seekIntersection(RaySt, Ray, node->kid[1], ob_2); dis_2 = dis_2 > EPS ? dis_2 : DBL_MAX;
	ob = dis_1 < dis_2 ? ob_1 : ob_2;
	return std::min(dis_1, dis_2);
}

double ObjectTree::seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object& ob) {
	switch (ob.type) {
	case PLANE:		return RayPlane		(RaySt, Ray,(*(Mat<>*)ob.v[0])[0], (*(Mat<>*)ob.v[0])[1], (*(Mat<>*)ob.v[0])[2], *(double*)ob.v[1]);
	case CIRCLE:	return RayCircle	(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[2], *(Mat<>*)ob.v[1]);
	case TRIANGLE:	return RayTriangle	(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>* )ob.v[1], *(Mat<>*)ob.v[2]);
	case PLANESHAPE:return RayPlaneShape(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>* )ob.v[1], *(Mat<>*)ob.v[2], (bool(*)(double, double))ob.v[3]);
	case SPHERE:	return RaySphere	(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[1], (bool(*)(double, double))ob.v[2]);
	case CUBOID:	return RayCuboid	(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>* )ob.v[1]);
	}
}
/*--------------------------------[ add Object ]--------------------------------*/
void ObjectTree::addPlane(Mat<>& n, Mat<>& p0, Material* material) {
	Object ob; ob.type = PLANE;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = n;(*(Mat<>*)ob.v[0]).normalize();
	ob.v[1] = new double;	*(double*)ob.v[1] = -((*(Mat<>*)ob.v[0])[0] * p0[0] + (*(Mat<>*)ob.v[0])[1] * p0[1] + (*(Mat<>*)ob.v[0])[2] * p0[2]);
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addPlane(std::initializer_list<double> nt, std::initializer_list<double> p0t, Material* material) {
	Mat<> n(3), p0(3); 
	n = nt;
	p0 = p0t;

	Object ob; ob.type = PLANE;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = n; (*(Mat<>*)ob.v[0]).normalize();
	ob.v[1] = new double;	*(double*)ob.v[1] = -((*(Mat<>*)ob.v[0])[0] * p0[0] + (*(Mat<>*)ob.v[0])[1] * p0[1] + (*(Mat<>*)ob.v[0])[2] * p0[2]);
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addCircle(Mat<>& center, double R, Mat<>& n, Material* material) {
	Object ob; ob.type = CIRCLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new Mat<>;	*(Mat<>*) ob.v[1] = n; (*(Mat<>*)ob.v[1]).normalize();
	ob.v[2] = new double;	*(double*)ob.v[2] = R;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addCircle(std::initializer_list<double> centert, double R, std::initializer_list<double> nt, Material* material) {
	Mat<> n(3), center(3);
	n = nt;
	center = centert;

	Object ob; ob.type = CIRCLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new Mat<>;	*(Mat<>*) ob.v[1] = n; (*(Mat<>*)ob.v[1]).normalize();
	ob.v[2] = new double;	*(double*)ob.v[2] = R;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addTriangle(Mat<>& p1, Mat<>& p2, Mat<>& p3, Material* material) {
	Object ob; ob.type = TRIANGLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = p1;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = p2;
	ob.v[2] = new Mat<>;	*(Mat<>*)ob.v[2] = p3;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addTriangle(std::initializer_list<double> p1t, std::initializer_list<double> p2t, std::initializer_list<double> p3t, Material* material) {
	Mat<> p1(3), p2(3), p3;
	p1 = p1t;
	p2 = p2t;
	p3 = p3t;

	Object ob; ob.type = TRIANGLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = p1;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = p2;
	ob.v[2] = new Mat<>;	*(Mat<>*)ob.v[2] = p3;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addPlaneShape(Mat<>& n, Mat<>& p0, bool(*f)(double, double), Material* material) {
	Object ob; ob.type = PLANESHAPE;
	ob.v = (void**)calloc(4, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = p0;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = n;  (*(Mat<>*)ob.v[1]).normalize();
	ob.v[2] = new Mat<>;	
	ob.v[3] = (void*)f;
	{
		if (n[0] == 0 && n[1] == 0)*(Mat<>*)ob.v[2] = { 1,0,0 };
		else {
			Mat<> t(3);
			(*(Mat<>*)ob.v[2]).cross_(*(Mat<>*)ob.v[1], t = { 0,0,1 }).normalize();
		}
	}
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addPlaneShape(std::initializer_list<double> nt, std::initializer_list<double> p0t, bool(*f)(double, double), Material* material) {
	Mat<> n(3), p0(3); 
	n = nt;
	p0 = p0t;

	Object ob; ob.type = PLANESHAPE;
	ob.v = (void**)calloc(4, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = p0;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = n;  (*(Mat<>*)ob.v[1]).normalize();
	ob.v[2] = new Mat<>;
	ob.v[3] = (void*)f;
	{
		if (n[0] == 0 && n[1] == 0)*(Mat<>*)ob.v[2] = { 1,0,0 };
		else {
			Mat<> t(3);
			(*(Mat<>*)ob.v[2]).cross_(*(Mat<>*)ob.v[1], t = { 0,0,1 }).normalize();
		}
	}
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addSphere(Mat<>& center, double r, Material* material, bool(*f)(double, double)) {
	Object ob; ob.type = SPHERE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new double;	*(double*)ob.v[1] = r;
	ob.v[2] = (void*)f;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addSphere(std::initializer_list<double> centert, double r, Material* material, bool(*f)(double, double)) {
	Mat<> center(3);
	center = centert;

	Object ob; ob.type = SPHERE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new double;	*(double*)ob.v[1] = r;
	ob.v[2] = (void*)f;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addCuboid(Mat<>& pmin, Mat<>& pmax, Material* material){
	Object ob; ob.type = CUBOID;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = pmin;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = pmax;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addCuboid(std::initializer_list<double> pmint, std::initializer_list<double> pmaxt, Material* material) {
	Mat<> pmin(3), pmax(3);
	pmin = pmint;
	pmax = pmaxt;

	Object ob; ob.type = CUBOID;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = pmin;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = pmax;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addStl(const char* file, Mat<>& center, double size, Material** material) {
	Mat<> p0(3), p1(3), p2(3), p3(3), p4(3), p5(3), p6(3); Mat<short> a;
	GraphicsFileCode::stlRead(file, p0, p1, p2, p3, a);
	for (int i = 0; i < p0.cols; i++) {
		addTriangle(
			((p4 = { p1(0,i), p1(1,i), p1(2,i) }) *= size) += center,
			((p5 = { p2(0,i), p2(1,i), p2(2,i) }) *= size) += center,
			((p6 = { p3(0,i), p3(1,i), p3(2,i) }) *= size) += center,
			material[a[i]]
		);
	}
}

void ObjectTree::addStl(const char* file, std::initializer_list<double> centert, double size, Material** material) {
	Mat<> center(3);
	center = centert;

	Mat<> p0(3), p1(3), p2(3), p3(3), p4(3), p5(3), p6(3); Mat<short> a;
	GraphicsFileCode::stlRead(file, p0, p1, p2, p3, a);
	for (int i = 0; i < p0.cols; i++) {
		addTriangle(
			((p4 = { p1(0,i), p1(1,i), p1(2,i) }) *= size) += center,
			((p5 = { p2(0,i), p2(1,i), p2(2,i) }) *= size) += center,
			((p6 = { p3(0,i), p3(1,i), p3(2,i) }) *= size) += center,
			material[a[i]]
		);
	}
}

/*#############################################################################

*						光线追踪  Ray Tracing

##############################################################################*/
/*--------------------------------[ 初始化 ]--------------------------------*/
void RayTracing::init(int width, int height) {
	ScreenPix.zero(height, width);
	Screen.   zero(height, width);
	for (int i = 0; i < Screen.size(); i++) Screen[i].zero(3);
}
/*--------------------------------[ 画像素 ]--------------------------------*/
void RayTracing::setPix(int x, int y, Mat<>& color) {
	if (x < 0 || x >= ScreenPix.rows || y < 0 || y >= ScreenPix.cols) return;
	ScreenPix(ScreenPix.rows - x - 1, y).R = std::min((int)(color[0] * 0xFF), 0xFF);
	ScreenPix(ScreenPix.rows - x - 1, y).G = std::min((int)(color[1] * 0xFF), 0xFF);
	ScreenPix(ScreenPix.rows - x - 1, y).B = std::min((int)(color[2] * 0xFF), 0xFF);
}
/*--------------------------------[ 渲染 ]--------------------------------
*	[过程]:
		[1] 计算屏幕矢量、屏幕X,Y向轴
		[2] 对屏幕每个像素遍历
			[3] 计算像素矢量、光线矢量、光线追踪起点
			[4] 光线追踪算法
			[5] 基于结果绘制该像素色彩
-------------------------------------------------------------------------*/
void RayTracing::paint(const char* fileName, int sampleSt, int sampleEd) {
	//[0]
	objTree.build();
	//[1]
	static Mat<> ScreenVec, ScreenXVec, ScreenYVec(3);
	ScreenVec. sub(gCenter, Eye);															//屏幕轴由眼指向屏幕中心
	ScreenYVec = { ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0], 1, 0 };			//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.cross(ScreenVec, ScreenYVec.normalize()).normalize();						//屏幕X向轴与屏幕轴、屏幕Y向轴正交

	//[2]
	static Mat<> PixYVec, PixXVec, PixVec, Ray, RaySt, color(3); clock_t start;
	for (int sample = sampleSt; sample < sampleEd; sample++) {
		double rate = 1.0 / (sample + 1), randX = RAND_DBL, randY = RAND_DBL;

		Screen *= 1 - rate;

		for (int x = 0; x < Screen.rows; x++) {
			for (int y = 0; y < Screen.cols; y++) {
				PixVec.add(														//[3]
					PixXVec.mul(x + randX - Screen.rows / 2 - 0.5, ScreenXVec),
					PixYVec.mul(y + randY - Screen.cols / 2 - 0.5, ScreenYVec)
				); 
				traceRay(														//[4][5]
					RaySt.add(gCenter,   PixVec), 
					Ray.  add(ScreenVec, PixVec).normalize(), 
					color.zero(), 0
				); 
				Screen(x, y) += (color *= rate);
			} 
		} if (sample % 1 == 0) printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));


		if (sample % 1 == 0) {
			for (int x = 0; x < Screen.rows; x++) 
				for (int y = 0; y < Screen.cols; y++) 
					setPix(x, y, Screen(x, y));

			GraphicsFileCode::ppmWrite(fileName, ScreenPix); start = clock();
		}
	}
}
/******************************************************************************
*						追踪光线
*	[步骤]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
&	[注]:distance > 1而不是> 0，是因为反射光线在接触面的精度内，来回碰自己....
******************************************************************************/
Mat<>& RayTracing::traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level) {
	//搜索光线-对象的最近交点距离
	Object* obj; 
	double dis = objTree.seekIntersection(RaySt, Ray, obj); 

	//无交点
	if (dis == DBL_MAX)
		color.zero();

	//光源
	else if (obj->material->rediate) 
		color = obj->material->color;

	//最大追踪层数
	else if (level > maxRayLevel) 
		return color;

	else {
		Material* material = obj->material;

		//计算面矢、交点
		static Mat<> faceVec(3), tmp, tmp2;

		RaySt += (tmp.mul(dis, Ray));

		switch (obj->type) {
		case PLANE:		faceVec = *(Mat<>*)obj->v[0]; break;
		case CIRCLE:	faceVec = *(Mat<>*)obj->v[1]; break;
		case TRIANGLE:	
			faceVec.cross_(
				tmp .sub(*(Mat<>*)obj->v[1], *(Mat<>*)obj->v[0]),
				tmp2.sub(*(Mat<>*)obj->v[2], *(Mat<>*)obj->v[0])
			).normalize(); break;
		case PLANESHAPE:faceVec = *(Mat<>*)obj->v[1]; break;
		case SPHERE:	faceVec.sub(RaySt, *(Mat<>*)obj->v[0]).normalize(); break;
		case CUBOID:	 
			if      (fabs(RaySt[0] - (*(Mat<>*)obj->v[0])[0]) < EPS || fabs(RaySt[0] - (*(Mat<>*)obj->v[1])[0]) < EPS) faceVec = { 1, 0, 0 };
			else if (fabs(RaySt[1] - (*(Mat<>*)obj->v[0])[1]) < EPS || fabs(RaySt[1] - (*(Mat<>*)obj->v[1])[1]) < EPS) faceVec = { 0, 1, 0 };
			else if (fabs(RaySt[2] - (*(Mat<>*)obj->v[0])[2]) < EPS || fabs(RaySt[2] - (*(Mat<>*)obj->v[1])[2]) < EPS) faceVec = { 0, 0, 1 };
			break;
		}

		//色散补丁
		static int refractColorIndex; 
		static double refractRateBuf; 
		static bool isChromaticDisperson;
		if (level == 0) 
			refractColorIndex = RAND_DBL * 3, refractRateBuf = 1, isChromaticDisperson = 0;
		
		static Mat<> Ray0; 
		Ray0 = Ray;

		//快速反射
		if (material->quickReflect) {								//Reflect Quick: 计算点光源直接照射该点产生的颜色
			double lightCos = 0, t;
			faceVec *= faceVec.dot(Ray) > 0 ? -1 : 1;
			for (int i = 0; i < PointLight.size(); i++) {
				t = faceVec.dot(tmp.sub(PointLight[i], RaySt).normalize());
				lightCos = t > lightCos ? t : lightCos;
			} 
			color = 1;
			color *= lightCos;
		}

		//漫反射
		else if (material->diffuseReflect) {
			diffuseReflect(Ray0, faceVec, Ray);
			traceRay(RaySt, Ray, color, level + 1);
			color *= material->reflectLossRate;
		}

		//反射
		else if (RAND_DBL < material->reflect) {
			reflect(Ray0, faceVec, Ray);
			traceRay(RaySt, Ray, color, level + 1);
			color *= material->reflectLossRate;
		}

		//折射
		else {
			//折射率补丁
			if (material->refractRate[0] != material->refractRate[1] 
			 || material->refractRate[0] != material->refractRate[2]) 
				isChromaticDisperson = 1;
			double t = refractRateBuf; 
			refractRateBuf = refractRateBuf == material->refractRate[refractColorIndex] ? 
				1 : material->refractRate[refractColorIndex];

			refract(Ray0, faceVec, Ray, t, refractRateBuf);
			traceRay(RaySt += (tmp.mul(EPS, Ray)), Ray, color, level + 1);
			color *= material->refractLossRate;
		}
		//色散补丁
		if (level == 0 && isChromaticDisperson) { 
			double t = color[refractColorIndex]; 
			color.zero();
			color[refractColorIndex] = 3 * t;
		}

		color.elementMul(obj->material->color);
	}

	//雾
	if(haze) 
		Haze(color, color, haze_A, dis, haze_beta);

	return color;
}
