#ifndef _PLOT_H
#define _PLOT_H
#include "Graphics.h"
#include <stdlib.h>
#include <stdio.h>
#include "Mat.h"
class Plot
{
public:
	Graphics* g = NULL;														//核心图形学类
	double pSizeMax[2], pSizeMin[2];									//XY轴范围
	double pDelta[2];													//单位一的像素格数
	/*---------------- SET 设置 ----------------*/
	Plot(void) { init(); }
	void init();											//初始化
	void clear(RGB color);												//清屏
	void setAxisRange(const double minx, const double miny, const double maxx, const double maxy);//设置轴范围
	/*---------------- 实数坐标 To 像素坐标 ----------------*/
	int coor2pix(double coor, int dim);									//坐标To像素坐标
	double pix2coor(int pix, int dim);
	int value2pix(double value, int dim);								//值To像素值
	/*---------------- PLOT ----------------*/
	void plotPoint(const double x, const double y);						//画点
	void plotLine(const double x1, const double y1, const double x2, const double y2);//画线
	void plotWave(const double x[], const double y[], const int n);		//画曲线
	void plotCircle(const double x, const double y, const double r);	//画圆
	void plotEllipse(const double x, const double y, const double rx, const double ry);//画椭圆
	void plotRectangle(const double sx, const double sy, const double ex, const double ey);//画矩形
	void plotVector(const double sx, const double sy, const double ex, const double ey);//画箭头
	void contour(Mat<double>& map, const int N);					//画等高线
	void contourface(Mat<double>& map, const int N);				//画等高线2
	void grid();														//显示网格
	RGB colorlist(const int N, const int i, const int model);			//色谱
};
#endif