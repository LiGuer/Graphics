#include "Object.h"

using namespace Matrix;
using namespace Intersect;
using namespace ObjectLib;
/*#############################################################################
* 
*					对象/对象树
*
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
			sub((*(Mat<>*)bound->v[0]), *(Mat<>*)ob->v[0], p);
			add((*(Mat<>*)bound->v[1]), *(Mat<>*)ob->v[0], p);
			break;
		case TRIANGLE:
			for (int j = 0; j < 3; j++) {
				(*(Mat<>*)bound->v[0])[j] = std::min((*(Mat<>*)ob->v[0])[j], std::min((*(Mat<>*)ob->v[1])[j], (*(Mat<>*)ob->v[2])[j]));
				(*(Mat<>*)bound->v[1])[j] = std::max((*(Mat<>*)ob->v[0])[j], std::max((*(Mat<>*)ob->v[1])[j], (*(Mat<>*)ob->v[2])[j]));
			}
			break;
		case SPHERE:
			p = *(double*)ob->v[1];
			sub((*(Mat<>*)bound->v[0]), *(Mat<>*)ob->v[0], p);
			add((*(Mat<>*)bound->v[1]), *(Mat<>*)ob->v[0], p);
			break;
		case ELLIPSOID: 
			p = std::max(1.0 / (*(Mat<>*)bound->v[1])(0, 0), 1.0 / (*(Mat<>*)bound->v[1])(1, 1));
			p = sqrt(std::max(p(0), 1.0 / (*(Mat<>*)bound->v[1])(2, 2)));
			sub((*(Mat<>*)bound->v[0]), *(Mat<>*)ob->v[0], p);
			add((*(Mat<>*)bound->v[1]), *(Mat<>*)ob->v[0], p);
			break;
		case CUBOID: delete bound; bound = ob; break;
		case RING:
			p = *(double*)ob->v[1];
			sub((*(Mat<>*)bound->v[0]), *(Mat<>*)ob->v[0], p);
			add((*(Mat<>*)bound->v[1]), *(Mat<>*)ob->v[0], p);
			break;
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
		if ((*(Mat<>*)a.bound->v[0])[dim] != (*(Mat<>*)b.bound->v[0])[dim])
			return (*(Mat<>*)a.bound->v[0])[dim] < (*(Mat<>*)b.bound->v[0])[dim];
		return (*(Mat<>*)a.bound->v[1])[dim] < (*(Mat<>*)b.bound->v[1])[dim];
		});
	build(obSet, l, (l + r) / 2, node->kid[0]);
	build(obSet, (l + r) / 2 + 1, r, node->kid[1]);
}

/*--------------------------------[ 求交 ]--------------------------------*/
double ObjectTree::seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object*& ob) {
	double disMin = seekIntersection(RaySt, Ray, root, ob), dis_t;
	for (int i = 0; i < planeNum; i++) {
		dis_t = seekIntersection(RaySt, Ray, *ObNodeList[i].ob);
		if (dis_t > EPS && disMin > dis_t) { 
			ob = ObNodeList[i].ob; 
			disMin = dis_t;
		}
	}
	return disMin;
}

double ObjectTree::seekIntersection(Mat<>& RaySt, Mat<>& Ray, ObjectNode* node, Object*& ob) {
	if (node->ob != NULL) { 
		ob = node->ob; 
		return seekIntersection(RaySt, Ray, *node->ob); 
	}
	if (seekIntersection(RaySt, Ray, *node->bound) == DBL_MAX) 
		return DBL_MAX;

	Object* ob_1, * ob_2;
	double dis_1 = seekIntersection(RaySt, Ray, node->kid[0], ob_1); 
	dis_1 = dis_1 > EPS ? dis_1 : DBL_MAX;

	double dis_2 = seekIntersection(RaySt, Ray, node->kid[1], ob_2); 
	dis_2 = dis_2 > EPS ? dis_2 : DBL_MAX;

	ob = dis_1 < dis_2 ? ob_1 : ob_2;
	return std::min(dis_1, dis_2);
}

double ObjectTree::seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object& ob) {
	switch (ob.type) {
	case PLANE:
		return RayPlane(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[1]);
	case CIRCLE:
		return RayCircle(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[2], *(Mat<>*)ob.v[1]);
	case TRIANGLE:
		return RayTriangle(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>*)ob.v[1], *(Mat<>*)ob.v[2]);
	case PLANESHAPE:
		return RayPlaneShape(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>*)ob.v[1], *(Mat<>*)ob.v[2], (bool(*)(double, double))ob.v[3]);
	case SPHERE:
		return RaySphere(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[1], (bool(*)(double, double))ob.v[2]);
	case ELLIPSOID:
		return RayQuadric(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>*)ob.v[1]);
	case CUBOID:
		return RayCuboid(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>*)ob.v[1]);
	case RING:
		return inter(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[1], *(double*)ob.v[2]);
	}
}

/*--------------------------------[ add Object ]--------------------------------*/
void ObjectTree::addPlane(Mat<>& n, Mat<>& p0, Material* material) {
	Object ob; ob.type = PLANE;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = n; normalize((*(Mat<>*)ob.v[0]));
	ob.v[1] = new double;	*(double*)ob.v[1] = dot(*(Mat<>*)ob.v[0], p0);
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addPlane(std::initializer_list<double> nt, std::initializer_list<double> p0t, Material* material) {
	Mat<> n(3), p0(3);
	n = nt;
	p0 = p0t;

	Object ob; ob.type = PLANE;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = n; normalize(*(Mat<>*)ob.v[0]);
	ob.v[1] = new double;	*(double*)ob.v[1] = dot(*(Mat<>*)ob.v[0], p0);
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addCircle(Mat<>& center, double R, Mat<>& n, Material* material) {
	Object ob; ob.type = CIRCLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new Mat<>;	*(Mat<>*) ob.v[1] = n; normalize(*(Mat<>*)ob.v[1]);
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
	ob.v[1] = new Mat<>;	*(Mat<>*) ob.v[1] = n; normalize(*(Mat<>*)ob.v[1]);
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
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = n;  normalize(*(Mat<>*)ob.v[1]);
	ob.v[2] = new Mat<>;
	ob.v[3] = (void*)f;
	{
		if (n[0] == 0 && n[1] == 0)*(Mat<>*)ob.v[2] = { 1,0,0 };
		else {
			Mat<> t(3);
			cross_(*(Mat<>*)ob.v[2], *(Mat<>*)ob.v[1], t = { 0,0,1 });
			normalize(*(Mat<>*)ob.v[2]);
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
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = n;  normalize(*(Mat<>*)ob.v[1]);
	ob.v[2] = new Mat<>;
	ob.v[3] = (void*)f;
	{
		if (n[0] == 0 && n[1] == 0)*(Mat<>*)ob.v[2] = { 1,0,0 };
		else {
			Mat<> t(3);
			cross_(*(Mat<>*)ob.v[2], *(Mat<>*)ob.v[1], t = { 0,0,1 });
			normalize(*(Mat<>*)ob.v[2]);
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

void ObjectTree::addEllipsoid(Mat<>& center, Mat<>& PInv, Material* material, bool(*f)(double, double)) {
	Object ob; ob.type = ELLIPSOID;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new Mat<>;	*(Mat<>*) ob.v[1] = PInv;
	ob.v[2] = (void*)f;
	ob.material = material;
	ObjectSet.push_back(ob);
}

void ObjectTree::addCuboid(Mat<>& pmin, Mat<>& pmax, Material* material) {
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
	GraphicsIO::stlRead(file, p0, p1, p2, p3, a);
	for (int i = 0; i < p0.cols; i++) {
		addTriangle(
			add(p4, mul(p4, size, p4 = { p1(0,i), p1(1,i), p1(2,i) }), center),
			add(p5, mul(p5, size, p5 = { p2(0,i), p2(1,i), p2(2,i) }), center),
			add(p6, mul(p6, size, p6 = { p3(0,i), p3(1,i), p3(2,i) }), center),
			material[a[i]]
		);
	}
}

void ObjectTree::addStl(const char* file, std::initializer_list<double> centert, double size, Material** material) {
	Mat<> center(3);
	center = centert;

	Mat<> p0(3), p1(3), p2(3), p3(3), p4(3), p5(3), p6(3); Mat<short> a;
	GraphicsIO::stlRead(file, p0, p1, p2, p3, a);
	for (int i = 0; i < p0.cols; i++) {
		addTriangle(
			add(p4, mul(p4, size, p4 = { p1(0,i), p1(1,i), p1(2,i) }), center),
			add(p5, mul(p5, size, p5 = { p2(0,i), p2(1,i), p2(2,i) }), center),
			add(p6, mul(p6, size, p6 = { p3(0,i), p3(1,i), p3(2,i) }), center),
			material[a[i]]
		);
	}
}

void ObjectTree::addRing(Mat<>& center, double R, double r, Material* material) {
	Object ob; ob.type = RING;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new double;	*(double*)ob.v[1] = R;
	ob.v[2] = new double;	*(double*)ob.v[2] = r;
	ob.material = material;
	ObjectSet.push_back(ob);
}