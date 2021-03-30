#include "src/GraphicsND.h"
#include "LiGu_AlgorithmLib/Fractal.h"
#include <Windows.h>
int main() {
	GraphicsND g(1000, 1000);
	g.g.clear(0xFAE6CD);

	Mat<double> axis(3, 1), center(3, 1);
	{double t[] = { 500 ,500 ,0 }; center.getData(t); }
	g.translation(center);
	{double t[] = { 1 ,1 ,0 }; axis.getData(t); }
	{double t[] = { 0 ,0 ,0 }; center.getData(t); }
	g.rotate(axis, 3/3, center);
	g.g.PaintColor = 0; g.g.PaintSize = 0;
	g.drawAxis(100, 100, 100, 1);
	// Tree
	g.FaceColor = 0x21a366; g.g.PaintSize = 0;
	std::vector<Mat<double>> TreeSt, TreeEd;
	{
		Mat<double> tMat(3, 1);
		{double t[] = { 0,0,0 }; tMat.getData(t); } TreeSt.push_back(tMat);
		{double t[] = { 0,0,0 + 200 }; tMat.getData(t); } TreeEd.push_back(tMat);
	}
	Fractal::FractalTree3D(TreeSt, TreeEd, 6, (double)30 * 2 * PI / 360);
	for (int i = 0; i < TreeSt.size(); i++) {
		Mat<double> tmp;
		double Width = 0.1 * (tmp.add(TreeEd[i], TreeSt[i].negative(tmp))).norm();
		g.drawFrustum(TreeSt[i], TreeEd[i], Width, 0.7 * Width, 45,true,true);
	}
	g.FaceColor = 0xFF0011; g.g.PaintSize = 0;
	Mat<double> p1(3, 1), p2(3, 1), p3(3, 1);
	{double t[] = { 0,0,500 }; p1.getData(t); }
	{double t[] = { 0,-500,0 }; p2.getData(t); }
	{double t[] = { 0,600,0 }; p3.getData(t); } 
	g.drawTriangle(p1, p2, p3, true);

	g.FaceColor = 0x00FF11; g.g.PaintSize = 0;
	{double t[] = { 0,0,500 }; p1.getData(t); }
	{double t[] = { -500,0,0 }; p2.getData(t); }
	{double t[] = { 600,0,0 }; p3.getData(t); }
	g.drawTriangle(p1, p2, p3, true);
	g.g.writeImg("D:/LIGU.ppm");
}