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
	void drawPoint		(Mat<ARGB>& image, int x0, int y0);								//»­µã
	void drawLine		(Mat<ARGB>& image, int x1, int y1, int x2, int y2);				//»­Ïß
	void drawCircle		(Mat<ARGB>& image, int x0, int y0, int r);						//»­Ô²
	void drawEllipse	(Mat<ARGB>& image, int x0, int y0, int rx, int ry);				//»­ÍÖÔ²
	void drawRectangle	(Mat<ARGB>& image, int x1, int y1, int x2, int y2);				//»­¾ØĞÎ
	void drawPolygon	(Mat<ARGB>& image, int* x, int* y, int n);						//»­¶à±ßĞÎ
	void drawWave		(Mat<ARGB>& image, int* x, int* y, int n);						//»­ÇúÏß
	void drawBezier		(Mat<ARGB>& image, int* x, int* y, int n);						//»­±´Èû¶ûÇúÏß
	void drawGrid		(Mat<ARGB>& image, int sx, int sy, int ex, int ey, int dx, int dy);	//»­Íø¸ñ

	void fillRectangle	(Mat<ARGB>& image, int sx, int sy, int ex, int ey, ARGB color);	//Ìî³äµ¥É«
	void fillFlood		(Mat<ARGB>& image, int x0, int y0, ARGB color);					//·ºÀÄÌî³ä
	void fillPolygon	(Mat<ARGB>& image, int* x, int* y, int n);						//¶à±ßĞÎÌî³ä

	void drawChar		(Mat<ARGB>& image, int x0, int y0, char charac);				//ÏÔÊ¾×Ö·û
	void drawString		(Mat<ARGB>& image, int x0, int y0, const char* str);			//ÏÔÊ¾×Ö·û´®
	void drawNum		(Mat<ARGB>& image, int x0, int y0, fp64 num);					//ÏÔÊ¾Êı×Ö
}

#endif