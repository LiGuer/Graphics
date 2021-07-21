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
namespace GeometricalOptics {
	Mat<>& reflect			(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//反射
	Mat<>& refract			(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO);	//折射
	Mat<>& diffuseReflect	(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//漫反射
}
/*---------------- 光线追踪 ----------------*/
class RayTracing {
public:
	struct Material {															//材质
		Mat<> color{ 3 }, refractRate{ 3 };
		bool 
			rediate        = 0,
			quickReflect   = 0,
			diffuseReflect = 0;
		double
			reflect = 1, reflectLossRate = 1,
			refract = 0, refractLossRate = 1;
	};
	struct Triangle { Mat<> p[3]; Material* material = NULL; };					//三角形
	/*---------------- 基础参数 ----------------*/
	Mat<> Eye{ 3 }, gCenter{ 3 };
	Mat<RGB>	ScreenPix;
	Mat<Mat<>>	Screen;
	int maxRayLevel = 5;
	double ScreenXSize, ScreenYSize, eps = 1e-4;
	std::vector<Triangle> TriangleSet;											//三角形集
	std::vector<Mat<>>    PointLight;											//点光源集(QuickReflect专用)
	/*---------------- 底层 ----------------*/
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }					//构造函数
	void init (int width, int height);											//初始化
	void setPix(int x, int y, Mat<>& color);									//画像素
	/*---------------- DRAW ----------------*/
	void paint(const char* fileName, int sampleSt = 0, int sampleEd = 0x7FFFFFFF);		//渲染
	Mat<>& traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);					//追踪光线
	double seekIntersection				(Triangle& triangle, Mat<>& RaySt, Mat<>& Ray);	//求交点
	double seekIntersection_RaySphere	(Triangle& triangle, Mat<>& RaySt, Mat<>& Ray);	//求交点
	double seekIntersection_RayTriangle	(Triangle& triangle, Mat<>& RaySt, Mat<>& Ray);	//求交点
	//add
	void drawTriangle	(Mat<>& p1, Mat<>& p2, Mat<>& p3,	Material* material = NULL);				//画三角形
	void drawSphere		(Mat<>& center, double r,			Material* material = NULL);				//画球
};
#endif