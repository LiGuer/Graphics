#ifndef GRAPHICS_2D_H
#define GRAPHICS_2D_H

#include <algorithm>
#include <queue>
#include <vector>
#include "../../../Math/src/Matrix/Mat.h"
#include "./Geometry/bezier_curve.h"
#include "RGB.h"
#include "font.h"

using namespace std;

#define PI 3.141592653589

namespace Graphics {
	typedef long long	int64; 
	typedef float		fp32;
	typedef double		fp64;

	extern ARGB PaintColor;
	extern int  PaintSize, FontSize;

	/*-------------------------------- DRAW --------------------------------*/
	// Basic Geometry
	void drawPoint		(Mat<ARGB>& image, int x0, int y0);								//画点
	void drawLine		(Mat<ARGB>& image, int x1, int y1, int x2, int y2);				//画线
	void drawLine		(Mat<ARGB>& image, int* x, int* y, int n);						//画折线
	void drawCircle		(Mat<ARGB>& image, int x0, int y0, int r);						//画圆
	void drawEllipse	(Mat<ARGB>& image, int x0, int y0, int rx, int ry);				//画椭圆
	void drawTriangle	(Mat<ARGB>& image, int x1, int y1, 
										   int x2, int y2, int x3, int y3);				//画三角形
	void drawRectangle	(Mat<ARGB>& image, int x1, int y1, int x2, int y2);				//画矩形
	void fillRectangle	(Mat<ARGB>& image, int sx, int sy, int ex, int ey, ARGB color);	//填充矩形
	void drawRegularPolygon	
						(Mat<ARGB>& image, int  x, int  y, int l, int n, double a0 = 0);//画正多边形
	void drawPolygon	(Mat<ARGB>& image, int* x, int* y, int n);						//画多边形
	void fillPolygon	(Mat<ARGB>& image, int* x, int* y, int n);						//填充多边形
	void drawBezier		(Mat<ARGB>& image, vector<vector<double>>& points, int n);

	// Tessellation
	void drawGrid		(Mat<ARGB>& image, int sx, int sy, int ex, int ey, int dx, int dy);	//画网格

	// Fill
	void fillFlood		(Mat<ARGB>& image, int x0, int y0, ARGB color);					//泛滥填充

	// Text
	void drawChar		(Mat<ARGB>& image, int x0, int y0, char charac);				//显示字符
	void drawString		(Mat<ARGB>& image, int x0, int y0, const char* str);			//显示字符串
	void drawNum		(Mat<ARGB>& image, int x0, int y0, fp64 num);					//显示数字
}

#endif