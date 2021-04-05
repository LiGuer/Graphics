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
#define PI 3.141592653589

class RayTracing {
public:
	struct Material {															//材质
		unsigned int color = 0;
		double reflex, refraction;
	};
	struct Triangle {															//三角形
		Mat<double> p[3];
		Material* material = NULL;
	};
	/*---------------- 基础参数 ----------------*/
	Graphics g;																	//核心图形学类
	Mat<double> Eye, gCenter;
	int maxRayLevel = 10;
	std::vector<Triangle> TriangleSet;											//三角形集
	/*---------------- 底层 ----------------*/
	RayTracing() { ; }
	RayTracing(int width, int height, int Dim = 3) { init(width, height); }		//构造函数
	~RayTracing() { ; }															//析构函数
	void init(int width, int height);											//初始化
	/*---------------- DRAW ----------------*/
	void paint();																		//渲染
	unsigned int traceRay(Mat<double>& RaySt, Mat<double>& Ray, int level);
	double seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, double& RayFaceDistance, Mat<double>& intersection);	//求交点
};

#endif