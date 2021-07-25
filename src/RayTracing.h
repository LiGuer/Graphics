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
#define RAND_DBL (rand() / double(RAND_MAX))
/*---------------- 几何光学 ----------------*/
Mat<>& reflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//反射
Mat<>& refract(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO);	//折射
Mat<>& diffuseReflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);						//漫反射
/*---------------- 求交点 ----------------*/
double RayTriangle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与三角面
double RayPolygon	(Mat<>& RaySt, Mat<>& Ray, Mat<>* p,  int n);					//求交-射线与多边面
double RaySphere	(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R);			//求交-射线与球
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与长方体
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax);			//求交-射线与长方体 (轴对齐)
/*---------------- 光线追踪 ----------------*/
class RayTracing {
public:
	// 数据结构
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
	enum { Triangle = 0, Polygon, Sphere, Cuboid };
	struct Object { int type; Mat<>* p; double* v; Material* material = NULL; };//物体
	struct OctTree {															//八叉树
		std::vector<Object> ObjectSet;
		OctTree* kid[8], *parent = NULL;
		void buildTree() {

		}
	};
	//基础参数 
	Mat<> gCenter{ 3 }, Eye{ 3 };
	Mat<RGB>	ScreenPix;
	Mat<Mat<>>	Screen;
	int maxRayLevel = 5;
	double ScreenXSize, ScreenYSize, eps = 1e-4;
	std::vector<Object> ObjectSet;											//三角形集
	std::vector<Mat<>>    PointLight;											//点光源集(QuickReflect专用)
	//函数
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }					//构造函数
	void init (int width, int height);											//初始化
	void setPix(int x, int y, Mat<>& color);									//画像素
	void paint(const char* fileName, int sampleSt = 0, int sampleEd = 0x7FFFFFFF);		//渲染
	Mat<>& traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);					//追踪光线
	double seekIntersection (Object& ob, Mat<>& RaySt, Mat<>& Ray);				//求交点
	//add
	void addTriangle	(Mat<>& p1, Mat<>& p2, Mat<>& p3,	Material* material = NULL);	//+三角形
	void addSphere		(Mat<>& center, double r,			Material* material = NULL);	//+球
	void addCuboid		(Mat<>& pmin, Mat<>& pmax,			Material* material = NULL);	//+长方体
	void addStl(const char* file, Mat<>& center, double size, Material** material) {
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
};
#endif