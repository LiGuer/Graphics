#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <queue>
#include "font.h"
#include "Mat.h"
/******************************************************************************
*                    ������������޹ص���������
******************************************************************************/
typedef unsigned char  BOOLEAN;
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
/******************************************************************************
*                    Graphics �����ͼ��ѧ
******************************************************************************/
struct Graphics {
	/*---------------- �������� ----------------*/
	INT32U gSize[2] = { 2048,2048 };										//���ڳߴ�
	RGB* Map = NULL;														//ͼ(�ײ�)
	RGB PaintColor;															//������ɫ
	INT32U PaintSize = 0;													//���ʴ�С
	Mat<FP64> gM = { 3 };													//�任����
	/*---------------- ���� ----------------*/
	const RGB TRANSPARENT = 0xFFFFFFFF;										//RGB:͸��
	const INT32U FontSize = 16;												//�ַ���С
	/*---------------- �ײ� ----------------*/
	~Graphics() { free(Map); }												//��������
	void init();															//��ʼ��
	void clear();	 														//����
	void clear(RGB color);	 												//����
	void setPoint(INT32U x, INT32U y);										//�ײ㻭��
	void setPoint(INT32U x, INT32U y, RGB color);							//�ײ㻭��(ָ����ɫ)
	RGB  readPoint(INT32U x, INT32U y); 									//���� 
	void PicWrite(const CHAR* filename);									//��ͼ(�ײ�)
	/*---------------- DRAW ----------------*/
	void drawPoint(INT32U x0, INT32U y0);									//����
	void drawLine(INT32U x1, INT32U y1, INT32U x2, INT32U y2);				//����
	void drawCircle(INT32U x0, INT32U y0, INT32U r);					    //��Բ
	void drawEllipse(INT32U x0, INT32U y0, INT32U rx,INT32U ry);			//����Բ
	void drawTriangle(INT32U x1, INT32U y1, INT32U x2, INT32U y2, INT32U x3, INT32U y3);//��������
	void drawRectangle(INT32U x1, INT32U y1, INT32U x2, INT32U y2);		   	//������
	void drawWave(INT32U x[], INT32U y[], INT32U n);						//������
	void drawBezier(INT32S x[], INT32S y[], INT32U n);						//������������
	void drawGrid(INT32U sx, INT32U sy, INT32U ex, INT32U ey, INT32U dx, INT32U dy);
	void drawCopy(INT32U x0, INT32U y0, Graphics* gt);						//���Ʊ��ͼ
	void fill(INT32U sx, INT32U sy, INT32U ex, INT32U ey);		   			//��䵥ɫ
	void fill(INT32U sx, INT32U sy, INT32U ex, INT32U ey, RGB color);		//��䵥ɫ(ָ����ɫ)
	void floodfill(INT32U x0, INT32U y0, RGB color);						//�������
	void drawChar(INT32U x0, INT32U y0, CHAR charac);						//��ʾ�ַ�
	void drawString(INT32U x0, INT32U y0, const CHAR* str, INT32U n);		//��ʾ�ַ���
	void drawNum(INT32U x0, INT32U y0, FP64 num);							//��ʾ����
	/*---------------- ��ά�任 TRANSFORMATION ----------------*/
	void translation(INT32S dx, INT32S dy);									//ƽ��
	void rotate(FP64 theta);												//��ת
	void rotate(FP64 theta, INT32S x0, INT32S y0);							//��ת(���ڻ�׼��)
	void scaling(FP64 sx, FP64 sy);											//����(>1ֱ����ɱ任)
	void reflect();															//����
	void shear();															//����
	void confirmTrans();													//ȷ�ϱ任
	/*---------------- SET ----------------*/
	bool judgeOutRange(INT32U x0, INT32U y0);								//�ж������Ƿ����
	void setGSize(INT32U width, INT32U height);								//���ô��ڳߴ�
};
#endif
