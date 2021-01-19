#include"LiGu_Graphics/Graphics3D.h"
#include <stdio.h>

int main() {
	Graphics3D g(1000, 1000);

	Mat<double> axis(3, 1), center(3, 1);
	{double t[] = { 1 ,0 ,0 }; axis.getData(t); }
	{double t[] = { 500 ,1000 ,0 }; center.getData(t); }
	g.rotate(axis, 1, center);
	{double t[] = { 0 ,1 ,0 }; axis.getData(t); }
	g.rotate(axis, 0.4, center);


	Mat<double> z(40, 40);
	for (int y = 0; y < 40; y++) {
		for (int x = 0; x < 40; x++) {
			double x0 = (x - 20) / 1;
			double y0 = (y - 20) / 1;
			z(x, y) = 100 * pow(sqrt(fabs(x0)) + sqrt(fabs(y0)), 0.75) + 100;
		}
	}
	g.drawSurface(z, 100, 900, 100, 900);

	g.g->PaintColor = 0xFF0000;
	for (int y = 0; y < 40; y++) {
		for (int x = 0; x < 40; x++) {
			double x0 = ((double)x - 20) / 20;
			double y0 = ((double)y - 20) / 20;
			if (1 - y0 * y0 - x0 * x0 < 0) z(x, y) = 0;
			else z(x, y) = 100 * sqrt(1 - y0 * y0 - x0 * x0);
		}
	}
	g.drawSurface(z, 100, 900, 100, 900);


	g.g->PaintColor = 0x00FF00;
	for (int y = 0; y < 40; y++) {
		for (int x = 0; x < 40; x++) {
			double x0 = ((double)x - 20) / 7;
			double y0 = ((double)y - 20) / 7;
			z(x, y) = 30 * sin(y0 * y0 + x0 * x0) + 100;
		}
	}
	g.drawSurface(z, 600, 1000, 600, 1000);
	g.g->PicWrite("D:/LIGU.ppm");
}
