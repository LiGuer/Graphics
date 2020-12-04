#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <queue>
#include "font.h"
#include "../LiGu_AlgorithmLib/Mat.h"
/******************************************************************************
*                    ������������޹ص���������
******************************************************************************/
typedef bool  BOOL;
typedef unsigned char  INT8U;			/* Unsigned  8 bit quantity       */
typedef signed   char  INT8S;			/* Signed    8 bit quantity       */
typedef unsigned short INT16U;			/* Unsigned 16 bit quantity       */
typedef signed   short INT16S;			/* Signed   16 bit quantity       */
typedef unsigned int   INT32U;			/* Unsigned 32 bit quantity       */
typedef signed int     INT32S;			/* Signed   32 bit quantity       */
typedef long long      INT64S;			/* Signed   64 bit quantity       */
typedef float          FP32;			/* Single precision floating point*/
typedef double         FP64;			/* Double precision floating point*/
typedef char    CHAR;			/* �ַ�       */
typedef INT32U  RGB;
typedef INT8U  RGBBASIC;
/******************************************************************************
*                    Graphics �����ͼ��ѧ
******************************************************************************/
class Graphics {
public:
	/*---------------- �������� ----------------*/
	INT32S gWidth = 100, gHeight = 100;										//���ڳߴ�
	RGBBASIC* Map = NULL;														//ͼ
	RGB PaintColor = 0xFFFFFF;												//������ɫ
	INT32S PaintSize = 0, FontSize=16;										//���ʴ�С//�ַ���С
	Mat<FP64> gM;															//�任����
	/*---------------- ���� ----------------*/
	const RGB TRANSPARENT = 0xFFFFFFFF;										//RGB:͸��
	/*---------------- �ײ� ----------------*/
	Graphics() { ; }
	Graphics(INT32S width, INT32S height) { init(width, height); }
	~Graphics() { free(Map); }												//��������
	void init();															//��ʼ��
	void init(INT32S width, INT32S height);									//��ʼ��
	void clear(RGB color);	 												//����
	void setPoint(INT32S x, INT32S y, RGB color);							//�ײ㻭��
	RGB  readPoint(INT32S x, INT32S y); 									//���� 
	void PicWrite(const CHAR* filename);									//��ͼ
	void confirmTrans();													//ȷ�ϱ任
	/*---------------- DRAW ----------------*/
	void drawPoint(INT32S x0, INT32S y0);									//����
	void drawLine(INT32S x1, INT32S y1, INT32S x2, INT32S y2);				//����
	void drawCircle(INT32S x0, INT32S y0, INT32S r);					    //��Բ
	void drawEllipse(INT32S x0, INT32S y0, INT32S rx,INT32S ry);			//����Բ
	void drawTriangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2, INT32S x3, INT32S y3);//��������
	void drawRectangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		   	//������
	void drawWave(INT32S x[], INT32S y[], INT32S n);						//������
	void drawBezier(INT32S x[], INT32S y[], INT32S n);						//������������
	void drawGrid(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy);//������
	void drawCopy(INT32S x0, INT32S y0, RGBBASIC* gt, INT32S width, INT32S height);//���Ʊ��ͼ
	void fill(INT32S sx, INT32S sy, INT32S ex, INT32S ey, RGB color);		//��䵥ɫ
	void fillflood(INT32S x0, INT32S y0, RGB color);						//�������
	void fillPolygon(INT32S x[], INT32S y[], INT32S n);						//��������
	void drawChar(INT32S x0, INT32S y0, CHAR charac);						//��ʾ�ַ�
	void drawString(INT32S x0, INT32S y0, const CHAR* str, INT32S n);		//��ʾ�ַ���
	void drawNum(INT32S x0, INT32S y0, FP64 num);							//��ʾ����
	/*---------------- ��ά�任 TRANSFORMATION ----------------*/
	void translation(INT32S dx, INT32S dy);									//ƽ��
	void rotate(FP64 theta, INT32S x0, INT32S y0);							//��ת
	void scaling(FP64 sx, FP64 sy);											//����(>1ֱ����ɱ任)
	/*---------------- SET ----------------*/
	BOOL judgeOutRange(INT32S x0, INT32S y0);								//�ж������Ƿ����
	void setSize(INT32S width, INT32S height);								//���ô��ڳߴ�
};
#endif
