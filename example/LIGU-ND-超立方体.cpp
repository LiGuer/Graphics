#include "src/GraphicsND.h"
#include "LiGu_AlgorithmLib/Fractal.h"
#include "src/DigitalImageProcessing.h"
int main() {
	GraphicsND g(1000, 1000, 4);

	double alpha = 0, belta = 0;
	while (true) {
		g.g.writeImg("D:/LIGU.ppm");
		g.clear(0);
		//rotate
		Mat<double> axis(3, 1), center(4, 1);
		{double t[] = { alpha ,belta ,1 }; axis.getData(t); }
		{double t[] = { 500 ,500 ,500 ,500 }; center.getData(t); }
		g.rotate(axis, alpha, center); alpha += 0.02; belta += 0.004;
		//drawAxis
		g.g.PaintColor = 0; g.g.PaintSize = 0;
		g.drawAxis(100, 100, 100, 1);
		//drawRectangle
		g.g.PaintColor = 0x0; g.FaceColor = 0xFF1111;
		Mat<double> p1(4, 1), p2(4, 1), p3(4, 1), p4(4, 1), p5(4, 1);
		{double t[] = { 600 ,600 ,600,600 }; p1.getData(t); }
		{double t[] = { 400 ,400 ,400,400 }; p2.getData(t); }
		{double t[] = { 410 ,410 ,410,410 }; p3.getData(t); }
		{double t[] = { 350 ,350 ,350,350 }; p4.getData(t); }
		{double t[] = { 50 ,50 ,50,50 }; p5.getData(t); }
		//g.drawSuperCuboid(p1, p2 ,1 ,1);
		//g.drawSuperCuboid(p3, p4 ,1 ,1);
		g.g.PaintColor = 0xFFFFFF;
		//g.drawSuperSphere(center, 300);
		g.drawGrid(p5, p1, p2);
	}
}