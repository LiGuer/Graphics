#ifndef _PLOT_H
#define _PLOT_H
#include "Graphics.h"
#include <stdlib.h>
#include <stdio.h>
#include "Mat.h"
class Plot
{
public:
	Graphics* g;														//����ͼ��ѧ��
	double pSizeMax[2], pSizeMin[2];									//XY�᷶Χ
	double pDelta[2];													//��λһ�����ظ���
	/*---------------- SET ���� ----------------*/
	void init(Graphics* gt);											//��ʼ��
	void clear(RGB color);												//����
	void setAxisRange(const double minx, const double miny, const double maxx, const double maxy);//�����᷶Χ
	/*---------------- ʵ������ To �������� ----------------*/
	int coor2pix(double coor, int dim);									//����To��������
	int value2pix(double value, int dim);								//ֵTo����ֵ
	/*---------------- PLOT ----------------*/
	void plotPoint(const double x, const double y);						//����
	void plotLine(const double x1, const double y1, const double x2, const double y2);//����
	void plotWave(const double x[], const double y[], const int n);		//������
	void plotCircle(const double x, const double y, const double r);	//��Բ
	void plotEllipse(const double x, const double y, const double rx, const double ry);//����Բ
	void plotRectangle(const double sx, const double sy, const double ex, const double ey);//������
	void plotVector(const double sx, const double sy, const double ex, const double ey);//����ͷ
	void contour(const Mat<double> *map, const int N);					//���ȸ���
	void contourface(const Mat<double>* map, const int N);				//���ȸ���2
	void grid();														//��ʾ����
	RGB colorlist(const int N, const int i, const int model);			//ɫ��
};
#endif