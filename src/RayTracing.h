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
#include "RGB.h"
#include "GraphicsND.h"
#include <vector>
#define PI 3.141592653589
#define RAND_DBL (rand() / double(RAND_MAX))
class RayTracing {
public:
	struct Material {															//材质
		Mat<double> color{ 3, 1 };
		bool rediateRate = 0, diffuseReflect = 0, quickReflect = 0;
		double reflectRate = 1, refractRate = 0;
	};
	struct Triangle { Mat<double> p[3];	Material* material = NULL; };			//三角形
	/*---------------- 基础参数 ----------------*/
	Mat<double> Eye{ 3,1 }, gCenter{ 3,1 };
	Mat<RGB> ScreenPix;
	Mat<Mat<double>> Screen;
	int SamplesNum = 1e9, maxRayLevel = 5;
	double refractRateBuffer = 1, eps = 1e-4, maxRayLevelPR = 1.0 / 6;
	std::vector<Triangle> TriangleSet;											//三角形集
	std::vector<Material> MaterialSet;											//材质集
	std::vector<Mat<double>> PointLight;										//点光源集(QuickReflect专用)
	/*---------------- 底层 ----------------*/
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }					//构造函数
	~RayTracing() { ; }															//析构函数
	void init(int width, int height);											//初始化
	void readImg(const char* fileName);											//读图
	void writeImg(const char* filename);										//写图
	void setPix(int x, int y, Mat<double>& color);								//画像素
	void readObj(const char* filename, Mat<double>& origin, Material* material = NULL);											//读3D模型
	/*---------------- DRAW ----------------*/
	void paint(const char* fileName, int sampleSt = 0);							//渲染
	Mat<double>& traceRay(Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& color, int level);
	double seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray);													//求交点
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