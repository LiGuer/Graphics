#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include "../LiGu_Codes/LiGu_Graphics/src/Fractal.h"
#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsFileCode.h"
#include <ctime>
void drawSurface(Mat<> z, double xs, double xe, double ys, double ye, GraphicsND* G) {
	Mat<> p(3), pl(3), pu(3), FaceVec, tmp, light(3); light.fill(1).normalized();
	double dx = (xe - xs) / z.rows, 
		   dy = (ye - ys) / z.cols;
	for (int y = 0; y < z.cols; y++) {
		for (int x = 0; x < z.rows; x++) {
			p.getData(xs + x * dx, ys + y * dy, z(x, y));
			for (int k = 0; k < 2; k++) {
				if ((k == 1 && (x == 0 || y == 0)) 
				||  (k == 0 && (x == z.rows - 1 || y == z.cols - 1))) continue;
				int dt = k == 0 ? 1 : -1;
				pl.getData(xs + (x + dt) * dx, ys + y * dy, z(x + dt, y));
				pu.getData(xs + x * dx, ys + (y + dt) * dy, z(x, y + dt));
				FaceVec.crossProduct(
					FaceVec.sub(pl, p),
					tmp.    sub(pu, p)
				).normalized();
				double t = (FaceVec.dot(light) + 1) / 2;
				G->FaceColor = (int)(t * 0xFF) * 0x10000 + (int)(t * 0xFF) * 0x100 + (int)(t * 0xFF);
				G->drawTriangle(p, pl, pu, 1, 0);		
			}
		}
	}
}
#define N 2000
int main() {
	GraphicsND G(1000, 1000);
	//旋转视角
	Mat<> t1(3), t2(3), t3(3), map(N, N), t = map;
	G.rotate(t1.getData(0, 0, 1), PI / 8, t2.zero());
	G.rotate(t1.getData(1, 0, 0), -PI / 3, t2.zero());
	//山海
	//map.functionIndex([](double x, int i) { return sin(sqrt(pow(i % N - N/2, 2) + pow(i / N - N / 2, 2)) / (2 * 10 * PI)); });
	//drawSurface(map *= 200, -400, 400, -400, 400, &G);
	
	srand(time(NULL));
	for (int k = 2; k <= map.rows; k *= 2) {
		Fractal::PerlinNoise(t, k);
		map += (t *= 1.0 / k);
	}
	map.function([](double x) { return x < 0 ? 0 : x; });
	drawSurface(map *= 800, -400, 400, -400, 400, &G);
	//人像
	Mat<> faceVec, p1, p2, p3, t0(3), light(3); light.fill(1).normalized();
	GraphicsFileCode::StlFileRead("D:/wamp64/www/LiGu/database/Art/雕塑/David.stl", faceVec, p1, p2, p3);
	p1 *= 0.5;
	p2 *= 0.5;
	p3 *= 0.5;
	for (int i = 0; i < faceVec.cols; i++) {
		faceVec.getCol(i, t0);
		p1.getCol(i, t1);
		p2.getCol(i, t2);
		p3.getCol(i, t3); 
		double t = (t0.dot(light)/ t0.norm() + 1) / 2;
		G.FaceColor = (int)(t * 0xFF) * 0x10000 + (int)(t * 0xFF) * 0x100 + (int)(t * 0xFF);
		G.drawTriangle(t1, t2, t3, 1, 0);
	}


	GraphicsFileCode::StlFileRead("D:/wamp64/www/LiGu/database/Art/雕塑/Venus.stl", faceVec, p1, p2, p3);
	p1 *= 0.8;
	p2 *= 0.8;
	p3 *= 0.8;
	for (int i = 0; i < faceVec.cols; i++) {
		faceVec.getCol(i, t0);
		p1.getCol(i, t1); t1[0] += 120;
		p2.getCol(i, t2); t2[0] += 120;
		p3.getCol(i, t3); t3[0] += 120;
		double t = (t0.dot(light) / t0.norm() + 1) / 2;
		G.FaceColor = (int)(t * 0xFF) * 0x10000 + (int)(t * 0xFF) * 0x100 + (int)(t * 0xFF);
		G.drawTriangle(t1, t2, t3, 1, 0);
	}
	G.g.writeImg("D:/ligu.ppm");
}