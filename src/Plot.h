#ifndef _PLOT_H
#define _PLOT_H
#include "Graphics.h"
#include <stdlib.h>
#include <stdio.h>
class Plot
{
public:
	Graphics* g;
	double pMax[2], pMin[2];											//XY�᷶Χ
	int WindowSize[2];													//�������ظ���
	double deltaXY[2];													//��λһ�����ظ���

	Plot(Graphics* gt);
	/*---------------- SET ----------------*/
	void clear();														//����
	void setXYRange(const double x[], const double y[], const int n);	//�Զ�����XY�᷶Χ
	/*---------------- PLOT ----------------*/
	void plotPoint(const double x, const double y);		//����
	void plotWave(const double x[], const double y[], const int n);		//������
	void plotCircle(const double x, const double y, const double r);	//��Բ
	void plotEllipse(const double x, const double y, const double rx, const double ry);//����Բ
	void plotRectangle(const double sx, const double sy, const double ex, const double ey);//������
	void grid();														//��ʾ����
	/*---------------- COOR TO PIX ----------------*/
	int coor2pix(double coor, int dim);									//����To��������
	int value2pix(double value, int dim);								//ֵTo����ֵ
};
#endif