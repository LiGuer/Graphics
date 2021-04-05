#include"src/RayTracing.h"

int main() {
	RayTracing ray(1000, 1000);
	
	{ double t[] = { 0,0,-1000 }; ray.Eye.getData(t); }
	{ double t[] = { 0,0,0 }; ray.gCenter.getData(t); }

	{
		RayTracing::Triangle triangle;
		Mat<double> point(3, 1);
		{ double t[] = { 0,0,100 }; point.getData(t); } triangle.p[0] = point;
		{ double t[] = { 200,0,100 }; point.getData(t); } triangle.p[1] = point;
		{ double t[] = { 100,-20,100 }; point.getData(t); } triangle.p[2] = point;
		ray.TriangleSet.push_back(triangle);

		{ double t[] = { 400,200,200 }; point.getData(t); } triangle.p[0] = point;
		{ double t[] = { 200,0,200 }; point.getData(t); } triangle.p[1] = point;
		{ double t[] = { 0,400,200 }; point.getData(t); } triangle.p[2] = point;
		ray.TriangleSet.push_back(triangle);

		{ double t[] = { -800,-800,0 }; point.getData(t); } triangle.p[0] = point;
		{ double t[] = { 1600,0,2000 }; point.getData(t); } triangle.p[1] = point;
		{ double t[] = { 0,800,800 }; point.getData(t); } triangle.p[2] = point;
		ray.TriangleSet.push_back(triangle);


		RayTracing::Material material;
		material.color = 0xFF0000; ray.MaterialSet.push_back(material);
		material.color = 0x00DD00; ray.MaterialSet.push_back(material);
		material.color = 0x0; ray.MaterialSet.push_back(material);

		ray.TriangleSet[0].material = &ray.MaterialSet[0];
		ray.TriangleSet[1].material = &ray.MaterialSet[1];
		ray.TriangleSet[2].material = &ray.MaterialSet[2];
	}

	ray.paint();
	ray.g.writeImg("D:/LIGU.ppm");

}