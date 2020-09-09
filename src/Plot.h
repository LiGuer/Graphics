#ifndef _PLOT_H
#define _PLOT_H
#include "Graphics.h"
#include <stdlib.h>
#include <stdio.h>
class Plot
{
public:
	Graphics* g;
	double pMax[2], pMin[2];											//XY轴范围
	int WindowSize[2];													//长宽像素格数
	double deltaXY[2];													//单位一的像素格数

	Plot(Graphics* gt);
	/*---------------- SET ----------------*/
	void clear();														//清屏
	void setXYRange(const double x[], const double y[], const int n);	//自动设置XY轴范围
	/*---------------- PLOT ----------------*/
	void plotPoint(const double x, const double y);		//画点
	void plotWave(const double x[], const double y[], const int n);		//画曲线
	void plotCircle(const double x, const double y, const double r);	//画圆
	void plotEllipse(const double x, const double y, const double rx, const double ry);//画椭圆
	void plotRectangle(const double sx, const double sy, const double ex, const double ey);//画矩形
	void grid();														//显示网格
	/*---------------- COOR TO PIX ----------------*/
	int coor2pix(double coor, int dim);									//坐标To像素坐标
	int value2pix(double value, int dim);								//值To像素值
};
#endif