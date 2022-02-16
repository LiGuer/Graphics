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
#ifndef RAY_TRACING_H
#define RAY_TRACING_H
#include <time.h>
#include <vector>
#include <algorithm>
#include "RGB.h"
#include "GraphicsFileCode.h"
#define PI 3.141592653589
#define EPS 10e-4
#define RAND_DBL (rand() / double(RAND_MAX))
/*---------------- 几何光学 ----------------*/
Mat<>& reflect			(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//反射
Mat<>& refract			(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO);	//折射
Mat<>& diffuseReflect	(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//漫反射

double Haze(double I, double A, double dis, double beta);										// 雾
Mat<>& Haze(Mat<>& I, Mat<>& O, Mat<>& A, double dis, double beta);

/*---------------- 求交点 ----------------*/
double RayPlane		(Mat<>& RaySt, Mat<>& Ray, double& A, double& B, double& C, double& D);	//求交-射线与平面
double RayCircle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, double& R, Mat<>& normal);	//求交-射线与圆
double RayTriangle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与三角面
double RayPolygon	(Mat<>& RaySt, Mat<>& Ray, Mat<>* p,  int n);					//求交-射线与多边面
double RayPlaneShape(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, Mat<>& normal, Mat<>& one, bool(*f)(double, double));//求交-射线与平面图形
double RaySphere	(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R, bool(*f)(double, double) = NULL);			//求交-射线与球
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与长方体
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax);			//求交-射线与长方体 (轴对齐)

/*---------------- 对象/对象树 ----------------*/
struct Material {															//材质
	Mat<> color{ 3 }, refractRate{ 3 };
	bool
		rediate = 0,
		quickReflect = 0,
		diffuseReflect = 0;
	double
		reflect = 1, reflectLossRate = 1,
		refract = 0, refractLossRate = 1;
};

enum { PLANE = 0, CIRCLE, TRIANGLE, POLTGON, PLANESHAPE, SPHERE, CUBOID };
struct Object { int type; void** v; Material* material = NULL; };		//物体
struct ObjectNode { Object* ob = NULL, * bound = NULL; ObjectNode* kid[2] = { NULL, NULL }; };
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
	void addPlane		(Mat<>& n, Mat<>& p0,				Material* material = NULL);	//+平面
	void addCircle		(Mat<>& center, double R, Mat<>& n,	Material* material = NULL);	//+圆
	void addTriangle	(Mat<>& p1,Mat<>& p2, Mat<>& p3,	Material* material = NULL);	//+三角形
	void addPlaneShape	(Mat<>& n, Mat<>& p0, bool(*f)(double,double), Material* material = NULL);	//+平面图形
	void addSphere		(Mat<>& center, double r,			Material* material = NULL, bool(*f)(double, double) = NULL);	//+球
	void addCuboid		(Mat<>& pmin, Mat<>& pmax,			Material* material = NULL);	//+长方体
	void addStl	(const char* file, Mat<>& center, double size, Material** material);


	void addPlane(std::initializer_list<double> n, std::initializer_list<double> p0, Material* material);
	void addCircle(std::initializer_list<double> center, double R, std::initializer_list<double> n, Material* material);
	void addTriangle(std::initializer_list<double> p1, std::initializer_list<double> p2, std::initializer_list<double> p3, Material* material);
	void addPlaneShape(std::initializer_list<double> n, std::initializer_list<double> p0, bool(*f)(double, double), Material* material);
	void addSphere(std::initializer_list<double> center, double r, Material* material, bool(*f)(double, double) = NULL);
	void addCuboid(std::initializer_list<double> pmin, std::initializer_list<double> pmax, Material* material);
	void addStl(const char* file, std::initializer_list<double> center, double size, Material** material);

};
/*---------------- 光线追踪 ----------------*/
class RayTracing {
public:
	//基础参数 
	Mat<> gCenter{ 3 }, Eye{ 3 };
	Mat<RGB>	ScreenPix;
	Mat<Mat<>>	Screen;
	int maxRayLevel = 6;
	double ScreenXSize, ScreenYSize;
	ObjectTree objTree;
	std::vector<Mat<>>  PointLight;													//点光源集(QuickReflect专用)

	bool haze = false; 
	Mat<> haze_A{ 3 };
	double haze_beta = 1;
	//函数
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }						//构造函数
	void init (int width, int height);												//初始化
	void setPix(int x, int y, Mat<>& color);										//画像素
	void paint(const char* fileName, int sampleSt = 0, int sampleEd = 0x7FFFFFFF);	//渲染
	Mat<>& traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);				//追踪光线
};

#endif