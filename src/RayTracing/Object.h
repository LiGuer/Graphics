#ifndef RAY_TRACING_OBJECT_H
#define RAY_TRACING_OBJECT_H

#include "Material.h"
#include "Intersect.h"

/*---------------- ����/������ ----------------*/
enum { PLANE = 0, CIRCLE, TRIANGLE, POLTGON, PLANESHAPE, SPHERE, CUBOID };

struct Object { 		//����
	int type; 
	void** v; 
	Material* material = NULL; 
};

struct ObjectNode {
	Object* ob = NULL, * bound = NULL; 
	ObjectNode* kid[2] = { NULL, NULL }; 
};

class  ObjectTree {
public:
	ObjectNode* root = NULL;
	ObjectNode* ObNodeList;
	std::vector<Object> ObjectSet;											//�����μ�
	int planeNum = 0;

	void build(std::vector<Object>& obSet);
	void build() { build(ObjectSet); };
	void build(ObjectNode* obSet, int l, int r, ObjectNode*& node);

	double seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object*& ob);
	double seekIntersection(Mat<>& RaySt, Mat<>& Ray, ObjectNode* node, Object*& ob);
	double seekIntersection(Mat<>& RaySt, Mat<>& Ray, Object& ob);

	//add
	void addPlane		(Mat<>& n, Mat<>& p0, Material* material = NULL);	//+ƽ��
	void addCircle		(Mat<>& center, double R, Mat<>& n, Material* material = NULL);	//+Բ
	void addTriangle	(Mat<>& p1,Mat<>& p2, Mat<>& p3, Material* material = NULL);	//+������
	void addPlaneShape	(Mat<>& n, Mat<>& p0, bool(*f)(double, double), Material* material = NULL);	//+ƽ��ͼ��
	void addSphere		(Mat<>& center, double r, Material* material = NULL, bool(*f)(double, double) = NULL);	//+��
	void addCuboid		(Mat<>& pmin, Mat<>& pmax, Material* material = NULL);	//+������
	void addStl			(const char* file, Mat<>& center, double size, Material** material);

	void addPlane		(std::initializer_list<double> n, std::initializer_list<double> p0, Material* material);
	void addCircle		(std::initializer_list<double> center, double R, std::initializer_list<double> n, Material* material);
	void addTriangle	(std::initializer_list<double> p1,std::initializer_list<double> p2, std::initializer_list<double> p3, Material* material);
	void addPlaneShape	(std::initializer_list<double> n, std::initializer_list<double> p0, bool(*f)(double, double), Material* material);
	void addSphere		(std::initializer_list<double> center, double r, Material* material, bool(*f)(double, double) = NULL);
	void addCuboid		(std::initializer_list<double> pmin, std::initializer_list<double> pmax, Material* material);
	void addStl			(const char* file, std::initializer_list<double> center, double size, Material** material);

};

#endif