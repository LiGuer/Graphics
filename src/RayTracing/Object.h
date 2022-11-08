#ifndef RAY_TRACING_OBJECT_H
#define RAY_TRACING_OBJECT_H

#include <vector>
#include <algorithm>
#include "Material.h"
#include "../../../LiGu_Math/src/Math/Geometry/Intersect.h"
#include "../GraphicsIO.h"
#include "C:\CppLiGu\QtLiGu\a.h"

#define EPS 10e-4

using namespace Matrix;

namespace ObjectLib {

/*---------------- 对象/对象树 ----------------*/
enum { PLANE = 0, CIRCLE, TRIANGLE, POLTGON, PLANESHAPE, SPHERE, CUBOID, ELLIPSOID, RING };

struct Object { 		//物体
	int type; 
	void** v; 
	Material* material = NULL; 
};

inline Mat<>& FaceVector(Object& obj, Mat<>& intersect, Mat<>& ans) {
	static Mat<> tmp1(3), tmp2(3);

	switch (obj.type) {
	case PLANE:		ans = *(Mat<>*)obj.v[0]; break;
	case CIRCLE:	ans = *(Mat<>*)obj.v[1]; break;
	case TRIANGLE:
		normalize(cross_(ans,
			sub(tmp1, *(Mat<>*)obj.v[1], *(Mat<>*)obj.v[0]),
			sub(tmp2, *(Mat<>*)obj.v[2], *(Mat<>*)obj.v[0])
		)); break;
	case PLANESHAPE:ans = *(Mat<>*)obj.v[1]; break;
	case SPHERE:	normalize(sub(ans, intersect, *(Mat<>*)obj.v[0])); break;
	case ELLIPSOID:
		ans = {
			intersect(0) * pow((*(Mat<>*)obj.v[1])(0, 0), 2),
			intersect(1) * pow((*(Mat<>*)obj.v[1])(1, 1), 2),
			intersect(2) * pow((*(Mat<>*)obj.v[1])(2, 2), 2)
		};
		normalize(ans);
		break;
	case CUBOID:
		if (fabs(intersect[0] - (*(Mat<>*)obj.v[0])[0]) < EPS 
		||  fabs(intersect[0] - (*(Mat<>*)obj.v[1])[0]) < EPS) 
			ans = { 1, 0, 0 };
		else if (fabs(intersect[1] - (*(Mat<>*)obj.v[0])[1]) < EPS 
			 ||  fabs(intersect[1] - (*(Mat<>*)obj.v[1])[1]) < EPS) 
			ans = { 0, 1, 0 };
		else if (fabs(intersect[2] - (*(Mat<>*)obj.v[0])[2]) < EPS 
			 ||  fabs(intersect[2] - (*(Mat<>*)obj.v[1])[2]) < EPS) 
			ans = { 0, 0, 1 };
		break;
	case RING:
		double 
			a = ((*(double*)obj.v[2]) * (*(double*)obj.v[2]) - (*(double*)obj.v[1]) * (*(double*)obj.v[1]))
			  + intersect[0] * intersect[0] + intersect[1] * intersect[1] + intersect[2] * intersect[2],
			b = 4 * (*(double*)obj.v[1]) * (*(double*)obj.v[1]);
		ans = {
			4 * intersect[0] * a - 2 * b * intersect[0],
			4 * intersect[1] * a - 2 * b * intersect[1],
			4 * intersect[2] * a,
		};
		normalize(ans);
		break;
	}
	return ans;
}

struct ObjectNode {
	Object* ob = NULL, * bound = NULL; 
	ObjectNode* kid[2] = { NULL, NULL }; 
};

class  ObjectTree {
public:
	ObjectNode* root = NULL;
	ObjectNode* ObNodeList;
	std::vector<Object> ObjectSet;											//三角形集
	int planeNum = 0;

	void build(std::vector<Object>& obSet);
	void build() { build(ObjectSet); };
	void build(ObjectNode* obSet, int l, int r, ObjectNode*& node);

	double seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object*& ob);
	double seekIntersection(Mat<>& RaySt, Mat<>& Ray, ObjectNode* node, Object*& ob);
	double seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object& ob);

	//add
	void addPlane		(Mat<>& n, Mat<>& p0, Material* material = NULL);	//+平面
	void addCircle		(Mat<>& center, double R, Mat<>& n, Material* material = NULL);	//+圆
	void addTriangle	(Mat<>& p1,Mat<>& p2, Mat<>& p3, Material* material = NULL);	//+三角形
	void addPlaneShape	(Mat<>& n, Mat<>& p0, bool(*f)(double, double), Material* material = NULL);	//+平面图形
	void addSphere		(Mat<>& center, double r, Material* material = NULL, bool(*f)(double, double) = NULL);	//+球
	void addEllipsoid	(Mat<>& center, Mat<>& PInv, Material* material = NULL, bool(*f)(double, double) = NULL);
	void addCuboid		(Mat<>& pmin, Mat<>& pmax, Material* material = NULL);	//+长方体
	void addStl			(const char* file, Mat<>& center, double size, Material** material);
	void addRing		(Mat<>& center, double R, double r, Material* material = NULL);

	void addPlane		(std::initializer_list<double> n, std::initializer_list<double> p0, Material* material);
	void addCircle		(std::initializer_list<double> center, double R, std::initializer_list<double> n, Material* material);
	void addTriangle	(std::initializer_list<double> p1,std::initializer_list<double> p2, std::initializer_list<double> p3, Material* material);
	void addPlaneShape	(std::initializer_list<double> n, std::initializer_list<double> p0, bool(*f)(double, double), Material* material);
	void addSphere		(std::initializer_list<double> center, double r, Material* material, bool(*f)(double, double) = NULL);
	void addCuboid		(std::initializer_list<double> pmin, std::initializer_list<double> pmax, Material* material);
	void addStl			(const char* file, std::initializer_list<double> center, double size, Material** material);

};

}

#endif