#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include <ctime>
void Delay(int time) { clock_t now = clock(); while (clock() - now < time); }//time*1000为秒数 
int main() {
	int Dim = 6;
	GraphicsND G(500, 500, Dim); G.g.PaintSize = 1;// G.LINE = 1; G.FACE = 0;
	Mat<> theta(Dim, Dim), p1(Dim), p2(Dim), zero(Dim);
	while (1) 
	{
		theta(0, 1) += 0.1;
		theta(0, 2) += 0.05;
		theta(0, 3) += 0.05;
		theta(0, 4) += 0.05;
		theta(0, 5) += 0.05;
		theta(1, 2) += 0.05;
		theta(1, 3) += 0.03;
		theta(2, 3) += 0.01;
		theta(2, 4) += 0.01;
		theta(2, 5) += 0.01;
		G.rotate(theta, zero); 
		G.g.PaintColor = 0xFFFFFF;
		G.drawSuperCuboid(p1 = 100, p2 = -100);
		G.g.writeImg("D:/LiGu.ppm"); G.clear(0); G.TransformMat.E(); Delay(100);
	}
}