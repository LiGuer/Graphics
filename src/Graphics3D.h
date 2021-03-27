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
#ifndef GRAPHICS3D_H
#define GRAPHICS3D_H
#include "Graphics.h"
const double PI = 3.141592653598;
class Graphics3D
{
public:
	/*---------------- 基础参数 ----------------*/
	Graphics g;																//核心图形学类
	Mat<int> Z_index;
	Mat<double> WindowSize{ 2,1 };											//窗口尺寸
	static Mat<double> TransformMat;										//变换矩阵
	/*---------------- 底层 ----------------*/
	Graphics3D() { ; }
	Graphics3D(int width , int height) { init(width, height); }				//构造函数
	~Graphics3D() { ; }														//析构函数
	void init(int width, int height);										//初始化
	void value2pix(int x0, int y0, int z0, int& x, int& y, int& z);			//  (2|3D)
	void value2pix(Mat<double>& p0, int& x, int& y, int& z);				//  (anyD)
	/*---------------- DRAW ----------------*/
	// 0-D
	void drawPoint(double x0, double y0, double z0 = 0);					//画点 (2|3D)
	void drawPoint(Mat<double>& p0);										//画点 (anyD)
	// 1-D
	void drawLine(double sx0, double ex0, double sy0, double ey0, double sz0 = 0, double ez0 = 0);//画直线 (2|3D)
	void drawLine(Mat<double>& sp, Mat<double>& ep);						//画直线 (anyD)
	void drawPolyline(Mat<double>* p, int n);								//画折线
	void contour(Mat<double>& map, const int N);							//画等高线
	// 2-D
	void drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3);			//画三角形
	void drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct = NULL);//画矩形
	void drawRectangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4);//画矩形
	void drawCircle(Mat<double>& center, double r, Mat<double>* direct = NULL);		//画圆
	void drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct = NULL);//画椭圆
	void drawPolygon(Mat<double> p[], int n);										//画多边形
	void drawSurface(Mat<double> z, double xs, double xe, double ys, double ye);	//画曲面
	void drawBezier();																//画贝塞尔曲面
	void contourface(Mat<double>& map, const int N);								//画等高面
	// 3-D
	void drawTetrahedron(Mat<double> p[]);									//画四面体
	void drawCuboid(Mat<double> pMin, Mat<double> pMax);					//画矩体
	void drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta = 5);//画圆台
	void drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta = 5);//画圆柱
	void drawSphere(Mat<double>& center, double r, int delta = 5);			//画球
	void drawSphere2(Mat<double>& center, double r, int n = 300);
	void drawEllipsoid(Mat<double>& center, Mat<double>& r);				//画椭球
	// Word
	void drawChar(Mat<double>& p0, char charac);							//显示字符
	void drawString(Mat<double>& p0, const char* str, int n);				//显示字符串
	void drawNum(Mat<double>& p0, double num);								//显示数字
	// Other
	void grid();															//显示网格
	Graphics::ARGB colorlist(double index, int model);						//色谱
	/*---------------- Transformation ----------------*/
	void translation(Mat<double>& delta, Mat<double>& transMat = TransformMat);	//平移
	void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& transMat = TransformMat);//三维旋转-四元数
	void scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat = TransformMat);//缩放
	/*---------------- SET ----------------*/
	void setAxisLim(Mat<double> pMin, Mat<double> pMax);					//判断坐标范围
};

#endif