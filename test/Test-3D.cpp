#include "LiGu_Graphics/Graphics3D.h"
#include <Windows.h>
int main() {
	Graphics3D g(1000, 1000);
	Mat<double> Axis(3, 1);
	Axis[0] = 0.7;
	Axis[1] = 0.4;
	Axis[2] = 0.59;

	Mat<double> center(3, 1);
	center[0] = 500;center[1] = 500;center[2] = 500;
	Mat<double> centerSphere(3, 1);
	centerSphere[0] = 400; centerSphere[1] = 500; centerSphere[2] = 300;

	Mat<double> points[4];
	for(int i=0;i<4;i++)points[i].rands(3, 1, 300, 900);

	double theta = 0;
	while(1){
		theta += 0.01;
		if (theta >= 2 * 3.1415926) theta = 0;

		g.TransformMat.E(4);
		g.rotate(Axis, theta, center);

		g.g->clear(0);
		g.drawSphere(centerSphere,200);
		g.g->PicWrite("D:/LIGU.ppm");
		Sleep(50);
	}
	/*
	while (true) {
		theta += 0.01;
		if (theta >= 2 * 3.1415926) theta = 0;
		
		g.TransformMat.E(4);
		g.rotate(Axis, theta, center);

		g.g->clear(0);
		g.drawTetrahedron(points);
		g.g->PicWrite("D:/LIGU.ppm");
		Sleep(100);
	}*/
	g.g->PicWrite("D:/LIGU.ppm");
}