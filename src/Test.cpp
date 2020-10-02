﻿#include "Graphics.h"
#include "Plot.h"
#include "math.h"
#include "Mat.h"

int main()
{
	Graphics* g = new Graphics;
	Plot* plot = new Plot;
	plot->init(g);

	const int N = 1000;
	double x[N], y[N];
	for (int i = 0; i < N; i++) {
		x[i] = (i - 500);
		y[i] = sin(x[i] * 2 * 3.14 / 360);
	}
	plot->setAxisRange(-500, -1, 500, 1);
	plot->grid();
	g->PaintColor = 0xFF0000;
	g->PaintSize = 2;
	plot->plotWave(x, y, N);
	g->PaintColor = 0xFFFF00;

	Mat<double> map;
	map.zero(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			map.data[i * map.cols + j] = i + j;
		}
	}
	g->PaintSize = 0;
	g->PaintColor = 0xFFFF00;
	plot->contourface(&map, 100);
	plot->contour(&map, 5);
	int x2[] = { 395,479,1199,1101 ,1294 ,1417,857 ,668 ,1111 };
	int y2[] = { 1887,1998,1433,1867, 1715 ,1171 ,1163 ,1314,1321 };
	g->fillPolygon(x2, y2, 9);

	g->PicWrite("D:/LIGU.ppm");
}

/***************************************************************************
int main()
{
	Graphics *g = new Graphics;
	g->setSize(2048, 2000);
	g->init();

	g->PaintColor = 0xFFFFFF;
	for (int i = 0; i < 50; i++) {
		g->PaintSize = i;
		g->drawPoint(50 + i * 38, 1800);
	}

	g->PaintSize = 5;
	g->PaintColor = 0xFFCC66;
	g->drawEllipse(1000, 1000, 500, 300);
	g->PaintColor = 0x0000FF;
	g->drawCircle(1000, 1000, 500);
	g->PaintSize = 3;
	g->PaintColor = 0x66FFCC;
	g->drawEllipse(1000, 1000, 200, 500);

	g->PaintColor = 0x6600CC;
	g->drawTriangle(100, 900, 300, 1500, 100, 1700);
	g->PaintColor = 0x00FFCC;
	g->drawRectangle(1900, 100, 1700, 300);
	g->fill(1850, 150, 1750, 250, 0x00FFCC);

	g->PaintColor = 0x0000FF;
	const int N = 100;
	int x[N], y[N];
	for (int i = 0; i < N; i++) {
		x[i] = i*10 + 500;
		y[i] = sin(0.02 * x[i])*100+200;
	}
	g->drawWave(x, y, N);

	g->PaintSize = 1;
	g->drawNum(1000, 1200, -1203.567);
	
	Graphics* gt = new Graphics;
	gt->setSize(500, 500);
	gt->init();
	gt->clear(gt->TRANSPARENT);
	gt->PaintSize = 2;
	gt->PaintColor = 0xFFCC66;
	gt->drawEllipse(100, 100, 50, 30);
	gt->PaintColor = 0x0000FF;
	gt->drawCircle(100, 100, 50);
	gt->PaintColor = 0x66FFCC;
	gt->drawEllipse(100, 100, 20, 50);
	gt->PaintColor = 0x00FFCC;
	gt->drawRectangle(190, 10, 160, 40);
	gt->fill(180, 20, 170, 30, 0x00FFCC);
	gt->rotate(3.14 * 60 / 360, 100, 100);
	gt->confirmTrans();
	g->drawCopy(10, 10, gt->Map, gt->gWidth, gt->gHeight);

	gt->PaintSize = 0;
	gt->PaintColor = 0xFF9900;
	gt->clear(gt->TRANSPARENT);
	gt->drawString(0, 0, "Ligu~", 5);
	gt->scaling(4, 3);
	g->drawCopy(0, 0, gt->Map, gt->gWidth, gt->gHeight);

	g->PaintColor = 0xFFEE00;
	g->PaintSize = 3;
	int x2[] = { 1900-200,1900,1900 +200,1900 ,1900 - 200 };
	int y2[] = { 600,600 - 500,600,600 - 100 ,600 };
	g->drawBezier(x2, y2, 5);
	g->floodfill(1900, 502, 0xFFEE00);

	g->PaintColor = 0x0055FF;
	g->PaintSize = 3;
	int x3[] = { 1800,1800 + 100,1800 + 250 ,1800 };
	int y3[] = { 1000,1000 - 150,1000  ,1000+ 230 };
	g->drawBezier(x3, y3, 4);
	int x4[] = { 1800 ,1800 - 100 ,1800 - 250  ,1800 };
	int y4[] = { 1000,1000 - 150,1000  ,1000 + 230 };
	g->drawBezier(x4, y4, 4);

	gt->PaintSize = 0;
	gt->PaintColor = 0x0055FF;
	gt->clear(gt->TRANSPARENT);
	gt->drawString(0, 0, "I Love U", 8);
	gt->scaling(2, 2.5);
	g->drawCopy(1800 - 120, 1000 + 100, gt->Map, gt->gWidth, gt->gHeight);


	gt->setSize(g->gWidth, g->gHeight);
	gt->init();
	gt->drawCopy(0, 0, g->Map, g->gWidth, g->gHeight);
	gt->scaling(0.2, 0.2);
	gt->confirmTrans();
	g->drawCopy(0, 300, gt->Map, gt->gWidth, gt->gHeight);

	g->fill(100, 100, 500, 500, 0x88FFFF00);
	g->fill(600, 600, 800, 800, 0x88FF00FF);
	g->fill(900, 900, 1200, 1200, 0x8800FFFF);
	g->fill(900, 1500, 1200, 1800, 0x8800FF00);
	g->fill(100, 900, 500, 1200, 0x880000FF);
	g->fill(600, 900, 800, 1200, 0x88FF0000);

	g->PicWrite("D:/LIGU.ppm");
}
***************************************************************************/