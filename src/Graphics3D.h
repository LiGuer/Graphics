#ifndef GRAPHICS3D_H
#define GRAPHICS3D_H
#include "Graphics.h"
struct GraphicsCoor {
	INT32S x, y, z;
};
class Graphics3D
{
public:
	/*---------------- �������� ----------------*/
	Graphics* g;															//����ͼ��ѧ��
	GraphicsCoor gSize;														//ͼ��С
	RGB PaintColor;															//������ɫ
	INT32S PaintSize = 0;													//���ʴ�С
	Mat<FP64> gM;															//�任����
	/*---------------- �ײ� ----------------*/
	void init();															//��ʼ��
	void clear();															//����
	/*---------------- DRAW ----------------*/
	void drawPoint(GraphicsCoor p0);										//����
	void drawLine(GraphicsCoor sp, GraphicsCoor ep);						//��ֱ��
	void drawPlane(GraphicsCoor p1, GraphicsCoor p2, GraphicsCoor p3);		//��ƽ��
	void drawSphere(GraphicsCoor p0, INT32S r);								//����
	void drawEllipsoid(GraphicsCoor p0, GraphicsCoor r);					//������
	void drawFace();														//������
	void drawBezier();														//������������
	void drawCopy(GraphicsCoor p0, RGB* gt, GraphicsCoor size);				//���Ʊ��3Dͼ
	void fill(GraphicsCoor sp, GraphicsCoor ep, RGB color);					//���
	void floodfill(GraphicsCoor p0, RGB color);								//�������
	void drawChar(GraphicsCoor p0, CHAR charac);							//��ʾ�ַ�
	void drawString(GraphicsCoor p0, const CHAR* str, INT32S n);			//��ʾ�ַ���
	void drawNum(GraphicsCoor p0, FP64 num);								//��ʾ����
	/*---------------- SET ----------------*/
	BOOL judgeOutRange(INT32S x0, INT32S y0);								//�ж������Ƿ����
	void setSize();															//���ô��ڳߴ�
};

#endif