#include"../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include <stdio.h>

int main() {
	GraphicsND g(1000, 1000);

	Mat<> axis(3), center(3);

	g.rotate(axis = { 1 ,0 ,0 }, 1, center = { 500 ,500 ,0 });
	g.rotate(axis = { 0 ,1 ,0 }, 0.4, center);
	
	Mat<> z(40, 40);
	for (int y = 0; y < 40; y++) {
		for (int x = 0; x < 40; x++) {
			double x0 = (x - 20) / 1;
			double y0 = (y - 20) / 1;
			z(x, y) = -400 * pow(sqrt(fabs(x0)) + sqrt(fabs(y0)), 0.75) + -400;
		}
	}
	g.drawSurface(z, -400, 400, -400, 400);

	g.g.PaintColor = 0xFF0000;
	for (int y = 0; y < 40; y++) {
		for (int x = 0; x < 40; x++) {
			double x0 = ((double)x - 20) / 20;
			double y0 = ((double)y - 20) / 20;
			if (1 - y0 * y0 - x0 * x0 < 0) z(x, y) = 0;
			else z(x, y) = -400 * sqrt(1 - y0 * y0 - x0 * x0);
		}
	}
	g.drawSurface(z, -400, 400, -400, 400);


	g.g.PaintColor = 0x00FF00;
	for (int y = 0; y < 40; y++) {
		for (int x = 0; x < 40; x++) {
			double x0 = ((double)x - 20) / 7;
			double y0 = ((double)y - 20) / 7;
			z(x, y) = 30 * sin(y0 * y0 + x0 * x0) + -400;
		}
	}
	g.drawSurface(z, 100, 500, 100, 500);
	g.g.writeImg("D:/LIGU.ppm");

}
