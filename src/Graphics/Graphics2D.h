#ifndef GRAPHICS_2D_H
#define GRAPHICS_2D_H

#include <algorithm>
#include <queue>
#include "../../../LiGu_Math/src/Math/Matrix/Mat.h"
#include "../RGB.h"
#include "../font.h"

namespace Graphics {
	typedef long long	int64; 
	typedef float		fp32;
	typedef double		fp64;

	extern ARGB PaintColor;
	extern int  PaintSize, FontSize;

	/*-------------------------------- DRAW --------------------------------*/
	void drawPoint		(Mat<ARGB>& image, int x0, int y0);								//����
	void drawLine		(Mat<ARGB>& image, int x1, int y1, int x2, int y2);				//����
	void drawCircle		(Mat<ARGB>& image, int x0, int y0, int r);						//��Բ
	void drawEllipse	(Mat<ARGB>& image, int x0, int y0, int rx, int ry);				//����Բ
	void drawRectangle	(Mat<ARGB>& image, int x1, int y1, int x2, int y2);				//������
	void drawPolygon	(Mat<ARGB>& image, int* x, int* y, int n);						//�������
	void drawWave		(Mat<ARGB>& image, int* x, int* y, int n);						//������
	void drawBezier		(Mat<ARGB>& image, int* x, int* y, int n);						//������������
	void drawGrid		(Mat<ARGB>& image, int sx, int sy, int ex, int ey, int dx, int dy);	//������

	void fillRectangle	(Mat<ARGB>& image, int sx, int sy, int ex, int ey, ARGB color);	//��䵥ɫ
	void fillFlood		(Mat<ARGB>& image, int x0, int y0, ARGB color);					//�������
	void fillPolygon	(Mat<ARGB>& image, int* x, int* y, int n);						//��������

	void drawChar		(Mat<ARGB>& image, int x0, int y0, char charac);				//��ʾ�ַ�
	void drawString		(Mat<ARGB>& image, int x0, int y0, const char* str);			//��ʾ�ַ���
	void drawNum		(Mat<ARGB>& image, int x0, int y0, fp64 num);					//��ʾ����
}

#endif