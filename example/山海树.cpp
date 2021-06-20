#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include "../LiGu_Codes/LiGu_Graphics/src/Fractal.h"
#include <ctime>

int main() {
	GraphicsND G(1000, 1000);
	//旋转视角
	Mat<> t1(3), t2(3), t3(3);
	G.rotate(t1 = { 0, 0, 1 }, PI / 8, t2.zero());
	G.rotate(t1 = { 1, 0, 0 }, -PI / 3, t2.zero());
	//箭头
	G.FaceColor = 0xFF0000;
	G.drawTriangle(
		t1 = { 0, 0, 50 },
		t2 = { 50, 0, 1 },
		t3 = { -50, 0, 1 }
	);
	G.FaceColor = 0x00FF00;
	G.drawTriangle(
		t1 = { 0, 0, 50 },
		t2 = { 0, 50, 1 },
		t3 = { 0, -50, 1 }
	);
	//山海
	srand(time(NULL));
	Mat<> map(400, 400), t = map;
	for (int k = 2; k <= map.rows; k *= 2)
		map += (Fractal::PerlinNoise(t, k) *= 1.0 / k);
	G.g.PaintColor = 0x00FF00; G.drawSurface(map *= 800, -500, 500, -500, 500);
	map.zero();
	G.g.PaintColor = 0x0066FF; G.drawSurface(map, -500, 500, -500, 500);
	//树
	std::vector<Mat<>> st, ed;
	st.push_back(t1 = { 0, 0, 0 });
	ed.push_back(t1 = { 0, 0, 150 });
	Fractal::FractalTree3D(st, ed, 6, (double)30 * 2 * PI / 360);
	G.FaceColor = 0xFFFF00; G.g.PaintColor = 0;
	for (int i = 0; i < st.size(); i++) {
		double Width = 0.1 * (t1.sub(st[i], ed[i])).norm();
		G.drawFrustum(st[i], ed[i], Width, 0.7 * Width, 45);
	}
	G.g.writeImg("D:/Ligu.ppm");
}