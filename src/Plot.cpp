#include "Plot.h"
/******************************************************************************
*                    SET 设置
******************************************************************************/
/*---------------- init ----------------*/
void Plot::init(Graphics* gt) {
	g = gt;
	g->init();
	pSizeMax[0] = pSizeMax[1] = pSizeMin[0] = pSizeMin[1] = 0;
}
/*---------------- clear ----------------*/
void Plot::clear(RGB color) {
	g->clear(color);
}
/*---------------- setAxisRange ----------------*/
void Plot::setAxisRange(const double minx, const double miny, const double maxx, const double maxy) {
	pSizeMax[0] = maxx; pSizeMax[1] = maxy;
	pSizeMin[0] = minx; pSizeMin[1] = miny;
	pDelta[0] = g->gWidth / (pSizeMax[0] - pSizeMin[0]);
	pDelta[1] = g->gHeight / (pSizeMax[1] - pSizeMin[1]);
}
/******************************************************************************
*                    实数坐标 To 像素坐标
******************************************************************************/
/*---------------- COOR TO PIX ----------------*/
int Plot::coor2pix(double coor, int dim) {
	return (coor - pSizeMin[dim]) * pDelta[dim];
}
/*---------------- value TO PIX ----------------*/
int Plot::value2pix(double value, int dim) {
	return value * pDelta[dim];
}
/******************************************************************************
*                    绘制 Plot
******************************************************************************/
/*---------------- 画点 ----------------*/
void Plot::plotPoint(const double x, const double y) {
	g->drawPoint(coor2pix(x, 0), coor2pix(y, 1));
}
/*---------------- 画曲线 ----------------*/
void Plot::plotWave(const double x[], const double y[], const int n) {
	INT32S* gx = (INT32S*)malloc(sizeof(INT32S) * n);
	INT32S* gy = (INT32S*)malloc(sizeof(INT32S) * n);
	for (int i = 0; i < n; i++) {
		gx[i] = coor2pix(x[i], 0);
		gy[i] = coor2pix(y[i], 1);
	}
	g->drawWave(gx, gy, n);
	free(gx);free(gy);
}
/*---------------- 画圆 ----------------*/
void Plot::plotCircle(const double x, const double y, const double r) {
	g->drawEllipse(coor2pix(x, 0), coor2pix(y, 1), value2pix(r, 0), value2pix(r, 1));
}
/*---------------- 画椭圆 ----------------*/
void Plot::plotEllipse(const double x, const double y, const double rx, const double ry) {
	g->drawEllipse(coor2pix(x, 0), coor2pix(y, 1), value2pix(rx, 0), value2pix(ry, 1));
}
/*---------------- 画矩形 ----------------*/
void Plot::plotRectangle(const double sx, const double sy, const double ex, const double ey) {
	g->drawRectangle(coor2pix(sx, 0), coor2pix(sy, 1), coor2pix(ex, 0), coor2pix(ey, 1));
}
/*---------------- 画等高线 ----------------*/
void Plot::contour(const Mat<double> *map, const int N)
{
	int x_step[] = { 1,0,1 }, y_step[] = { 0,1,1 };
	double max = map->max(), min = map->min();			//get the max & min of the map
	double delta = (max - min) / N, layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {		//for N layer between max & min, get the edge of each layer
		for (int y = 0; y < map->rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map->cols - 1; x++) {
				int flag = map->getValue(x, y) >= layer ? 1 : 0;
				for (int k = 0; k < 3; k++) {			//basic unit is 2x2 matrix
					int xt = x + x_step[k];
					int yt = y + y_step[k];
					int flagtemp = map->getValue(xt, yt) >= layer ? 1 : 0;
					if (flagtemp != flag) { flag = 2; break; }
				}
				if (flag == 2) {
					for (int k = 0; k < 3; k++) {
						int xt = x + x_step[k];
						int yt = y + y_step[k];
						if (map->getValue(xt, yt) >= layer) {
							g->drawPoint(xt, yt);
						}
					}
				}
			}
		}
	}
}
/*---------------- 画等高线2 ----------------*/
void Plot::contourface(const Mat<double>* map, const int N)
{
	int x_step[] = { 1,0,1 }, y_step[] = { 0,1,1 };
	double max = map->max(), min = map->min();			//get the max & min of the map
	double delta = (max - min) / N, layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {		//for N layer between max & min, get the edge of each layer
		for (int y = 0; y < map->rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map->cols - 1; x++) {
				if (map->data[x * map->cols + y] >= layer)
					g->setPoint(x, y, colorlist(N, i, 1));
			}
		}
	}
}
/*---------------- 画网格 ----------------*/
void Plot::grid() {
	int x0 = coor2pix(0, 0), y0 = coor2pix(0, 1);		//原点像素坐标
	double size[2] = { pSizeMax[0] - pSizeMin[0],pSizeMax[1] - pSizeMin[1] };
	/*------ 网格 ------*/
	g->PaintColor = 0xffcc00;
	g->PaintSize = 1;
	/*------ 计算间隔值 ------*/
	size[0] /= 10;size[1] /= 10;
	double delta[2] = { 1,1 };	
	for (int dim = 0; dim < 2; dim++) {
		while ((int)size[dim] == 0) {
			size[dim] *= 10;
			delta[dim] /= 10;
		}
		while ((int)size[dim] >= 10) {
			size[dim] /= 10;
			delta[dim] *= 10;
		}
		delta[dim] *= (int)size[dim];
	}
	g->drawGrid(x0, y0, 0, 0, -value2pix(delta[0], 0), -value2pix(delta[1], 1));
	g->drawGrid(x0, y0, g->gWidth, 0, value2pix(delta[0], 0), -value2pix(delta[1], 1));
	g->drawGrid(x0, y0, 0, g->gHeight, -value2pix(delta[0], 0), value2pix(delta[1], 1));
	g->drawGrid(x0, y0, g->gWidth, g->gHeight, value2pix(delta[0], 0), value2pix(delta[1], 1));
	/*------ 坐标轴 ------*/
	g->PaintColor = 0xffffff;
	g->PaintSize = 3;
	g->drawLine(x0, coor2pix(pSizeMin[1], 1), x0, coor2pix(pSizeMax[1], 1));
	g->drawLine(coor2pix(pSizeMin[0], 0), y0, coor2pix(pSizeMax[0], 0), y0);
	/*------ 轴标号 ------*/
	g->PaintSize = 1;
	g->FontSize = 50;
	g->PaintColor = 0xffffff;
	g->drawChar(x0 - 40, y0 + 15, '0');
}
/*---------------- 色谱 ----------------*/
RGB Plot::colorlist(const int N, const int i, const int model)
{
	double R = 0, G = 0, B = 0, a = (double)i / N, b = 1 - a;
	switch (model)
	{
	case 1: {
		B = a <= 9.0 / 16 ? (a < 1.0 / 16 ? 0.5 + 8 * a : (a > 6.0 / 16 ? 1 - (16/3.0) * (a - 6.0 / 16) : 1)) : 0;
		R = b <= 9.0 / 16 ? (b < 1.0 / 16 ? 0.5 + 8 * b : (b > 6.0 / 16 ? 1 - (16/3.0) * (b - 6.0 / 16) : 1)) : 0;
		G = (a >= 3.0 / 16 && b >= 3.0 / 16) ? (a < 6.0 / 16 ? (16 / 3.0) * (a - 3.0 / 16) : (b < 6.0 / 16 ? (16 / 3.0) * (b - 3.0 / 16) : 1)) : 0;
	}break;
	}
	R *= 0xFF, G *= 0xFF, B *= 0xFF;
	return (RGB)R + (RGB)G * 0x100 + (RGB)B * 0x10000;
}
