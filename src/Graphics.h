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
#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <queue>
#include "font.h"
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "RGB.h"
/******************************************************************************
*                    Graphics �����ͼ��ѧ
******************************************************************************/
class Graphics {
public:
	/*-------------------------------- �������޹��������� --------------------------------*/
	typedef signed   char  INT8S;			/* Signed    8 bit quantity       */
	typedef unsigned short INT16U;			/* Unsigned 16 bit quantity       */
	typedef signed   short INT16S;			/* Signed   16 bit quantity       */
	typedef unsigned int   INT32U;			/* Unsigned 32 bit quantity       */
	typedef signed int     INT32S;			/* Signed   32 bit quantity       */
	typedef long long      INT64S;			/* Signed   64 bit quantity       */
	typedef float          FP32;			/* Single precision floating point*/
	#define TRANSPARENT 0xFFFFFFFF
	/*-------------------------------- �������� --------------------------------*/
	Mat<RGB> Canvas{ 100, 100 };											//ͼ
	Mat<FP64> TransMat;														//�任����
	ARGB PaintColor = 0xFFFFFF;												//������ɫ
	INT32S PaintSize = 0, FontSize=16;										//���ʴ�С//�ַ���С
	/*-------------------------------- �ײ㺯�� --------------------------------*/
	Graphics() { ; }
	Graphics(INT32S width, INT32S height) { init(width, height); }
	~Graphics() { }															//��������
	void init(INT32S width = 100, INT32S height = 100);						//��ʼ��
	void clear(ARGB color);	 												//����
	void setPoint(INT32S x, INT32S y, ARGB color);							//�ײ㻭��
	ARGB readPoint(INT32S x, INT32S y); 									//���� 
	void readImg(const char* filename);										//��ͼ
	void writeImg(const char* filename);									//��ͼ
	bool judgeOutRange(INT32S x0, INT32S y0);								//�жϹ���
	void transSelf();														//ȫͼ�任
	void CutSelf(INT32S sx, INT32S sy, INT32S ex, INT32S ey);				//����ͼ
	/*-------------------------------- DRAW --------------------------------*/
	void drawPoint(INT32S x0, INT32S y0);									//����
	void drawLine(INT32S x1, INT32S y1, INT32S x2, INT32S y2);				//����
	void drawCircle(INT32S x0, INT32S y0, INT32S r);					    //��Բ
	void drawEllipse(INT32S x0, INT32S y0, INT32S rx,INT32S ry);			//����Բ
	void drawRectangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		   	//������
	void drawPolygon(INT32S x[], INT32S y[], INT32S n);						//�������
	void drawWave(INT32S x[], INT32S y[], INT32S n);						//������
	void drawBezier(INT32S x[], INT32S y[], INT32S n);						//������������
	void drawGrid(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy);//������
	void drawCopy(INT32S x0, INT32S y0, Mat<RGB>& gt);						//���Ʊ��ͼ
	void fillRectangle(INT32S sx, INT32S sy, INT32S ex, INT32S ey, ARGB color);	//��䵥ɫ
	void fillFlood(INT32S x0, INT32S y0, ARGB color);						//�������
	void fillPolygon(INT32S x[], INT32S y[], INT32S n);						//��������
	void drawChar(INT32S x0, INT32S y0, char charac);						//��ʾ�ַ�
	void drawString(INT32S x0, INT32S y0, const char* str, INT32S n);		//��ʾ�ַ���
	void drawNum(INT32S x0, INT32S y0, FP64 num);							//��ʾ����
	/*-------------------------------- ��ά�任 --------------------------------*/
	void translation(INT32S dx, INT32S dy);									//ƽ��
	void rotate(FP64 theta, INT32S x0, INT32S y0);							//��ת
	void scaling(FP64 sx, FP64 sy, INT32S x0, INT32S y0);					//����
};
#endif
