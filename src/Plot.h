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
#ifndef _PLOT_H
#define _PLOT_H
#include "Graphics.h"
#include <stdlib.h>
#include <stdio.h>
class Plot
{
public:
	Graphics* g = NULL;													//����ͼ��ѧ��
	double pSizeMax[2], pSizeMin[2];									//XY�᷶Χ
	double pDelta[2];													//��λһ�����ظ���
	/*---------------- SET ���� ----------------*/
	Plot(void) { init(); }
	void init();														//��ʼ��
	void clear(RGB color);												//����
	void setAxisRange(const double minx, const double miny, const double maxx, const double maxy);//�����᷶Χ
	/*---------------- ʵ������ To �������� ----------------*/
	int coor2pix(double coor, int dim);									//����To��������
	double pix2coor(int pix, int dim);
	int value2pix(double value, int dim);								//ֵTo����ֵ
	/*---------------- PLOT ----------------*/
	void plotPoint(const double x, const double y);						//����
	void plotLine(const double x1, const double y1, const double x2, const double y2);//����
	void plotWave(const double x[], const double y[], const int n);		//������
	void plotCircle(const double x, const double y, const double r);	//��Բ
	void plotEllipse(const double x, const double y, const double rx, const double ry);//����Բ
	void plotRectangle(const double sx, const double sy, const double ex, const double ey);//������
	void plotVector(const double sx, const double sy, const double ex, const double ey);//����ͷ
	void contour(Mat<double>& map, const int N);						//���ȸ���
	void contourface(Mat<double>& map, const int N);					//���ȸ���2
	void grid();														//��ʾ����
	RGB colorlist(const int N, const int i, const int model);			//ɫ��
};
#endif