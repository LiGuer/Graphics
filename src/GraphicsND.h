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
#ifndef GRAPHICS_ND_H
#define GRAPHICS_ND_H
#include "Graphics.h"
#define PI 3.141592653589
class GraphicsND
{
public:
	/*---------------- 基础参数 ----------------*/
	Graphics g;																//核心图形学类
	Mat<Mat<int>> Z_Buffer;
	Mat<double> WindowSize{ 2,1 };											//窗口尺寸
	static Mat<double> TransformMat;										//变换矩阵
	unsigned int FaceColor;
	/*---------------- 底层 ----------------*/
	GraphicsND() { ; }
	GraphicsND(int width , int height, int Dim = 3) { init(width, height , Dim); }//构造函数
	~GraphicsND() { ; }														//析构函数
	void init(int width, int height, int Dim = 3);							//初始化
	void clear(ARGB color);										//清屏
	void value2pix(int x0, int y0, int z0, int& x, int& y, int& z);			//点To像素 (<=3D)
	void value2pix(Mat<double>& p0, Mat<int>& pAns);						//点To像素 (anyD)
	bool setPix(int x, int y, int z = 0, int size = -1);					//写像素 (正投影) (<=3D)
	bool setPix(Mat<int>& p0, int size = -1);								//写像素 (正投影) (anyD)
	void setAxisLim(Mat<double>& pMin, Mat<double>& pMax);					//设置坐标范围
	/*---------------- DRAW ----------------*/
	// 0-D
	void drawPoint(double x0 = 0, double y0 = 0, double z0 = 0);			//画点 (<=3D)
	void drawPoint(Mat<double>& p0);										//画点 (anyD)
	// 1-D
	void drawLine(double sx0 = 0, double ex0 = 0, double sy0 = 0, double ey0 = 0, double sz0 = 0, double ez0 = 0);					//画直线 (<=3D)
	void drawLine(Mat<double>& sp0, Mat<double>& ep0);						//画直线 (anyD)
	void drawPolyline(Mat<double>* p, int n, bool close = false);			//画折线
	void drawPolyline(Mat<double>& y, double xmin, double xmax);			//画折线
	void drawBezierLine(Mat<double> p[], int n);							//画贝塞尔曲线
	// 2-D
	void drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, bool FACE = false, bool LINE = true);						//画三角形
	void drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct = NULL, bool FACE = false, bool LINE = true);			//画矩形
	void drawQuadrilateral(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, bool FACE = false, bool LINE = true);//画四边形
	void drawPolygon(Mat<double> p[], int n, bool FACE = false, bool LINE = true);													//画多边形
	void drawCircle(Mat<double>& center, double r, Mat<double>* direct = NULL, bool FACE = false, bool LINE = true);				//画圆
	void drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct = NULL);										//画椭圆
	void drawSurface(Mat<double> z, double xs, double xe, double ys, double ye);													//画曲面
	void drawBezierFace(Mat<double> p[], int n);																					//画贝塞尔曲面
	// 3-D
	void drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, bool FACE = false, bool LINE = true);	//画四面体
	void drawCuboid(Mat<double>& pMin, Mat<double>& pMax, bool FACE = false, bool LINE = true);										//画矩体
	void drawPolyhedron(Mat<double>* p, int n, bool FACE = false, bool LINE = true);												//画多面体
	void drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta = 5, bool FACE = false, bool LINE = true);//画圆台
	void drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta = 5, bool FACE = false, bool LINE = true);			//画圆柱
	void drawSphere(Mat<double>& center, double r, int delta = 5, bool FACE = false, bool LINE = true);								//画球
	void drawSphere2(Mat<double>& center, double r, int n = 300);																	//画球
	void drawEllipsoid(Mat<double>& center, Mat<double>& r);																		//画椭球
	void drawBody(Mat<double>& center, Mat<double>& r);																				//画曲体
	void drawBezierBody(Mat<double> p[], int n);																					//画贝塞尔曲体
	// Word
	void drawChar(Mat<double>& p0, char charac);							//显示字符
	void drawString(Mat<double>& p0, const char* str, int n);				//显示字符串
	void drawNum(Mat<double>& p0, double num);								//显示数字
	// any-D
	void drawSuperLine(Mat<double>* p0, bool FACE = false, bool LINE = true);														//画线 any-D
	void drawSuperCuboid(Mat<double>& pMin, Mat<double>& pMax, bool FACE = false, bool LINE = true);								//画立方体 any-D
	void drawSuperSphere(Mat<double>& center, double r, bool FACE = false, bool LINE = true);										//画球体 any-D
	void drawGrid(Mat<double>& delta, Mat<double>& max, Mat<double>& min, bool LINE = true);										//画网格
	// Other
	void drawAxis(double Xmax = 0, double Ymax = 0, double Zmax = 0, bool negative = false);										//画坐标轴
	void contour(Mat<double>& map, const int N);																					//画等高线
	void contourface(Mat<double>& map, const int N);																				//画等高面
	ARGB colorlist(double index, int model);																				//色谱
	/*---------------- Transformation ----------------*/
	static void translation(Mat<double>& delta, Mat<double>& transMat = TransformMat);														//平移
	static void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& transMat = TransformMat);					//旋转
	static void scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat = TransformMat);									//缩放
};
#endif