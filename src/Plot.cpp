#include "Plot.h"
Plot::Plot(Graphics* gt) {
	g = gt;
	WindowSize[0] = g->gSize[0];
	WindowSize[1] = g->gSize[1];
	pMax[0] = pMin[0] = pMax[1] = pMin[1] = 0;
}

void Plot::clear() {
	g->clear();
}

void Plot::setXYRange(const double x[], const double y[], const int n) {
	for (int i = 0; i < n; i++) {
		pMax[0] = pMax[0] > x[i] ? pMax[0] : x[i];
		pMin[0] = pMin[0] < x[i] ? pMin[0] : x[i];
		pMax[1] = pMax[1] > y[i] ? pMax[1] : y[i];
		pMin[1] = pMin[1] < y[i] ? pMin[1] : y[i];
	}
	deltaXY[0] = WindowSize[0] / (pMax[0] - pMin[0]);
	deltaXY[1] = WindowSize[1] / (pMax[1] - pMin[1]);
}

void Plot::plotPoint(const double x, const double y) {
	g->drawPoint(coor2pix(x, 0), coor2pix(y, 1));
}

void Plot::plotWave(const double x[], const double y[], const int n) {
	INT32S* gx = (INT32S*)malloc(sizeof(INT32S) * n);
	INT32S* gy = (INT32S*)malloc(sizeof(INT32S) * n);
	setXYRange(x, y, n);
	for (int i = 0; i < n; i++) {
		*(gx + i) = coor2pix(x[i], 0);
		*(gy + i) = coor2pix(y[i], 1);
	}
	g->drawWave(gx, gy, n);
	free(gx);
	free(gy);
}
void Plot::plotCircle(const double x, const double y, const double r) {
	int gx = coor2pix(x, 0);
	int gy = coor2pix(y, 1);
	int grx = value2pix(r, 0);
	int gry = value2pix(r, 1);
	g->drawEllipse(gx, gy, grx, gry);
}
void Plot::plotEllipse(const double x, const double y, const double rx, const double ry) {
	int gx = coor2pix(x, 0);
	int gy = coor2pix(y, 1);
	int grx = value2pix(rx, 0);
	int gry = value2pix(ry, 1);
	g->drawEllipse(gx, gy, grx, gry);
}
void Plot::plotRectangle(const double sx, const double sy, const double ex, const double ey) {
	int gsx = coor2pix(sx, 0);
	int gsy = coor2pix(sy, 1);
	int gex = coor2pix(ex, 0);
	int gey = coor2pix(ey, 1);
	g->drawRectangle(gsx, gsy, gex, gey);
}
void Plot::contour(const Mat<double> *map, const int N)
{
	double max = map->max(), min = map->min();		//get the max & min of the map
	double delta = (max - min) / N;
	int x_step[] = { 1,0,1 }, y_step[] = { 0,1,1 };
	double layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {//for N layer between max & min, get the edge of each layer
		for (int y = 0; y < map->rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map->cols - 1; x++) {
				int flag = map->getValue(x, y) >= layer ? 1 : 0;
				for (int k = 0; k < 3; k++) {		//basic unit is 2x2 matrix
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
RGB Plot::colorlist(const int N, const int i, const int model)
{
	RGB color = 0;
	double R = 0, G = 0, B = 0;
	double a = (double)i / N, b = 1 - a;
	switch (model)
	{
	case 1: {
		B = a <= 9.0 / 16 ? (a < 1.0 / 16 ? 0.5 + 8 * a : (a > 6.0 / 16 ? 1 - (16/3.0) * (a - 6.0 / 16) : 1)) : 0;
		R = b <= 9.0 / 16 ? (b < 1.0 / 16 ? 0.5 + 8 * b : (b > 6.0 / 16 ? 1 - (16/3.0) * (b - 6.0 / 16) : 1)) : 0;
		G = (a >= 3.0 / 16 && b >= 3.0 / 16) ? (a < 6.0 / 16 ? (16 / 3.0) * (a - 3.0 / 16) : (b < 6.0 / 16 ? (16 / 3.0) * (b - 3.0 / 16) : 1)) : 0;
	}break;
	case 2: {
		B = a <= 9.0 / 16 ? (a < 1.0 / 16 ? 0.5 + 8 * a : (a > 6.0 / 16 ? 1 - (16 / 3.0) * (a - 6.0 / 16) : 1)) : 0;
		R = b <= 9.0 / 16 ? (b < 1.0 / 16 ? 0.5 + 8 * b : (b > 6.0 / 16 ? 1 - (16 / 3.0) * (b - 6.0 / 16) : 1)) : 0;
		G = (a >= 3.0 / 16 && b >= 3.0 / 16) ? (a < 6.0 / 16 ? (16 / 3.0) * (a - 3.0 / 16) : (b < 6.0 / 16 ? (16 / 3.0) * (b - 3.0 / 16) : 1)) : 0;
	}break;
	}
	printf("%.3f - %.3f - %.3f\n", R, G, B);
	R *= 0xFF, G *= 0xFF, B *= 0xFF;
	return color = (RGB)R + (RGB)G * 0x100 + (RGB)B * 0x10000;
}
void Plot::contourface(const Mat<double>* map, const int N)
{
	RGB color;
	double max = map->max(), min = map->min();		//get the max & min of the map
	double delta = (max - min) / N;
	int x_step[] = { 1,0,1 }, y_step[] = { 0,1,1 };
	double layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {//for N layer between max & min, get the edge of each layer
		color = colorlist(N, i, 2);
		for (int y = 0; y < map->rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map->cols - 1; x++) {
				if (map->getValue(x, y) >= layer)
					g->setPoint(x, y, color);
			}
		}
	}
}
int Plot::coor2pix(double coor, int dim) {
	return (coor - pMin[dim]) * deltaXY[dim];
}

int Plot::value2pix(double value, int dim) {
	return value * deltaXY[dim];
}

void Plot::grid() {
	double size[2] = { pMax[0] - pMin[0],pMax[1] - pMin[1] };
	size[0] /= 10;
	size[1] /= 10; 

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
	}
	int gs[2];
	for (int dim = 0; dim < 2; dim++) {
		double s = (int)(pMin[dim] / delta[dim]) * delta[dim];
		gs[dim] = coor2pix(s, dim);
	}
	int ge[2];
	for (int dim = 0; dim < 2; dim++) {
		double e = (int)(pMax[dim] / delta[dim]) * delta[dim] + 1;
		ge[dim] = coor2pix(e, dim);
	}
	g->drawGrid(gs[0], gs[1], ge[0], ge[1], value2pix(delta[0],0), value2pix(delta[1],1));
}