#include "LiGu_Graphics/Graphics3D.h"

int main() {
	Graphics3D g(1000, 1000);
	Mat<double> Axis(3, 1);
	Axis[0] = 0.7;
	Axis[1] = 0.4;
	Axis[2] = 0.59;
	Mat<double> center(3, 1);
	center[0] = 500;
	center[1] = 500;
	center[2] = 500;
	Mat<double> points(3, 4);
	points.rands(3, 100, 300, 900);
	double theta = 0;
	while (true) {
		theta += 0.1;
		if (theta >= 2 * 3.1415926) break;

		g.TransformMat.E(4);
		g.rotate(Axis, theta, center);

		for (int i = 0; i < 4; i++) {
			Mat<double> temp(3, 1);
			Mat<double> temp2(3, 1);
			for (int j = 0; j < 3; j++)temp[j] = points(j, i);
			for (int j = 0; j < 3; j++)temp2[j] = points(j, (i+1)%3);
			g.drawLine(temp, temp2);
		}
	}
	g.g->PicWrite("D:/LIGU.ppm");
}