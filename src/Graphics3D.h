﻿/*
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
#ifndef GRAPHICS3D_H
#define GRAPHICS3D_H
#include "Graphics.h"
const double PI = 3.141592653598;
class Graphics3D
{
public:
	/*---------------- 基础参数 ----------------*/
	Graphics* g = NULL;														//核心图形学类
	Mat<FP64> WindowSize{ 2,1 };											//窗口尺寸
	Mat<FP64> TransformMat;													//变换矩阵
	/*---------------- 底层 ----------------*/
	Graphics3D(int WindowSize_Width, int WindowSize_height);				//构造函数
	~Graphics3D();															//析构函数
	void init(int WindowSize_Width, int WindowSize_height);					//初始化
	/*---------------- DRAW ----------------*/
	void drawPoint(Mat<double>& p0);										//画点
	void drawLine(Mat<double>& sp, Mat<double>& ep);						//画直线
	void drawPlane(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3);		//画平面
	void drawPolygon(Mat<double> p[], int n);								//画多边形
	void drawTetrahedron(Mat<double> p[]);									//画四面体
	void drawCuboid(Mat<double> pMin, Mat<double> pMax);					//画矩体
	void drawCircle(Mat<double>& center, double r, Mat<double>& direction);	//画圆
	void drawSphere(Mat<double>& center, double r);							//画球
	void drawEllipsoid(Mat<double>& center, Mat<double>& r);				//画椭球
	void drawFace();														//画曲面
	void drawBezier();														//画贝塞尔曲面
	void fill(Mat<double>& sp, Mat<double>& ep, RGB color);					//填充
	void fillTriangle(Mat<double> p0[]);									//三角填充
	void fillflood(Mat<double>& p0, RGB color);								//泛滥填充
	void drawChar(Mat<double>& p0, CHAR charac);							//显示字符
	void drawString(Mat<double>& p0, const CHAR* str, INT32S n);			//显示字符串
	void drawNum(Mat<double>& p0, FP64 num);								//显示数字
	/*---------------- Transformation ----------------*/
	void translation(Mat<double>& delta, Mat<double>& translationMat);		//平移
	void translation(Mat<double>& delta);									//平移
	void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& rotateMat);//三维旋转-四元数
	void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center);//三维旋转-四元数
	void scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& scaleMat);//缩放
	void scaling(Mat<double>& scale, Mat<double>& center);					//缩放
	/*---------------- SET ----------------*/
	bool judgeOutRange(INT32S x0, INT32S y0);								//判断坐标是否过界
	void setSize();															//设置窗口尺寸
};

#endif