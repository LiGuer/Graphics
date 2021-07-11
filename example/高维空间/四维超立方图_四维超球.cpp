#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include <ctime>
void Delay(int time) { clock_t now = clock(); while (clock() - now < time); }//time*1000为秒数 
int main() {
	GraphicsND G(1000, 1000, 4); G.g.PaintSize = 1;// G.LINE = 1; G.FACE = 0;
	Mat<> p1(4), p2(4), zero(4);
	double a = 0.01, b = 0.001;
	while (1) {
		a += 0.01; b += 0.001;
		G.rotate(p1, p2, a, b, zero); G.g.PaintColor = 0xFFFFFF;
		G.drawSuperCuboid(p1 = { 100 ,100 ,100,100 }, p2 = { -100 ,-100 ,-100,-100 });
		G.drawSuperSphere(zero, 300);
		//G.draw4DSphere(zero, 300);
		G.g.writeImg("D:/LiGu.ppm"); G.clear(0); G.TransformMat.E(); Delay(10);
	}
}