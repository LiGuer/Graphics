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

Reference.
[1]Introduction Algorithms.THOMAS H.CORMEN,CHARLES E.LEISERSON,RONALD L.RIVEST,CLIFFORD STEIN
==============================================================================*/
#ifndef COMPUTATIONAL_GEOMETRY_H
#define COMPUTATIONAL_GEOMETRY_H
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <vector>
#include <stack>
#include"../../LiGu_AlgorithmLib/Mat.h"
#define PI 3.141592653589
namespace Geometry {
	/*----------------[ 2D 二维 ]----------------*/
	//三角形
	bool inTriangle(Mat<double>& p0, Mat<double>& TriP1, Mat<double>& TriP2, Mat<double>& TriP3);//是否在三角内
	double RayTriIntersection(Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& TriP1, Mat<double>& TriP2, Mat<double>& TriP3);	//射线与三角形交点
	//圆
	bool onCircle(Mat<double> Points[]);													//判断四点共圆
	Mat<double>& ThreePoints2Circle(Mat<double> Points[], Mat<double>& center, double& R);	//平面三点确定圆方程
	//Other
	Mat<double>* ConvexHull(Mat<double> point[], int n, int& ansPointNum);					//凸包
	Mat<double>* Delaunay(Mat<double> point[], int n, int& TrianglesNum);					//Delaunay三角剖分
	/*----------------[ 3D 三维 ]----------------*/
	//球
	double RaySphereIntersection(Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& Center, double R);	//射线与球面交点
	Mat<double>* getSphereFibonacciPoint(int& n);											//球面均匀点分布
}
#endif