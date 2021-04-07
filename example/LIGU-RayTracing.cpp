#include"src/RayTracing.h"

int main() {
	RayTracing ray(1000, 1000);
	
	{ double t[] = { -2000,0,0 }; ray.Eye.getData(t); }
	{ double t[] = { 0,0,0 }; ray.gCenter.getData(t); }

	{
		RayTracing::Material* material;
		Mat<double> p1(3, 1), p2(3, 1), p3(3, 1), p4(3, 1), p5(3, 1);

		material = new RayTracing::Material;
		material->color = 0xFF0000;
		{ double t[] = { 100,-100,100 }; p1.getData(t); }
		{ double t[] = { 2000,100,200 }; p2.getData(t); } 
		ray.drawCuboid(p1,p2, material);


		material = new RayTracing::Material;
		material->color = 0x00FF00;
		{ double t[] = { 200,200,400 }; p1.getData(t); }
		{ double t[] = { 800,300,500 }; p2.getData(t); }
		ray.drawCuboid(p1, p2, material);

		material = new RayTracing::Material;
		material->color = 0xFFFFFF; material->reflexRate = 1;
		{ double t[] = { -2000,-2000,-500 }; p1.getData(t); }
		{ double t[] = { -2000,2000,-500 }; p2.getData(t); }
		{ double t[] = { 2000,2000,-500 }; p3.getData(t); }
		{ double t[] = { 2000,-2000,-500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);


		{ double t[] = { -1000,-1000,1000 }; p1.getData(t); }
		ray.LightSource.push_back(p1);
	}

	ray.paint();
	ray.g.writeImg("D:/LIGU.ppm");

}