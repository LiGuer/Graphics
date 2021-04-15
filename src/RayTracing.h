/*
Copyright 2020 LiGuer. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/
#ifndef RAY_TRACING_H
#define RAY_TRACING_H
#include "Graphics.h"
#include "GraphicsND.h"
#include "RGB.h"
#define PI 3.141592653589

class RayTracing {
public:
	struct Material {															//材质
		RGB color = 0;
		bool rediateRate = 0, diffuseReflect = 0;
		double reflectRate = 1, refractRate = 0;
	};
	struct Triangle { Mat<double> p[3];	Material* material = NULL; };			//三角形
	/*---------------- 基础参数 ----------------*/
	Graphics g;																	//核心图形学类
	Mat<double> Eye{ 3,1 }, gCenter{ 3,1 };
	int SamplesNum = 10, maxRayLevel = 6;
	double refractRateBuffer = 1, eps = 1e-4;
	std::vector<Triangle> TriangleSet;											//三角形集
	std::vector<Material> MaterialSet;											//材质集
	std::vector<Mat<double>> LightSource;										//材质集
	/*---------------- 底层 ----------------*/
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }					//构造函数
	~RayTracing() { ; }															//析构函数
	void init(int width, int height);											//初始化
	/*---------------- DRAW ----------------*/
	void paint();																//渲染
	RGB traceRay(Mat<double>& RaySt, Mat<double>& Ray, RGB& color, int level);
	double seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, Mat<double>& intersection);	//求交点
	//几何光学 Geometrical Optics
	static Mat<double>& reflect(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& reflectRay);								//反射
	static Mat<double>& refract(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& refractRay, double rateIn, double rateOut);//折射
	static Mat<double>& diffuseReflect(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& refractRay);						//漫反射
	//2-D
	void drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Material* material = NULL);							//画三角形
	void drawQuadrilateral(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, Material* material = NULL);		//画四边形
	void drawPolygon(Mat<double> p[], int n, Material* material = NULL);														//画多边形
	void drawSurface(Mat<double> z, double xs, double xe, double ys, double ye);												//画曲面	
	// 3-D
	void drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, Material* material = NULL);		//画四面体
	void drawCuboid(Mat<double>& pMin, Mat<double>& pMax, Material* material = NULL);											//画矩体
	void drawPolyhedron(Mat<double>* p, int n, Material* material = NULL);														//画多面体
	void drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta = 5, Material* material = NULL);	//画圆台
	void drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta = 5, Material* material = NULL);					//画圆柱
	void drawSphere(Mat<double>& center, double r, Material* material = NULL);													//画球
	void drawEllipsoid(Mat<double>& center, Mat<double>& r, Material* material = NULL);											//画椭球
};

#endif