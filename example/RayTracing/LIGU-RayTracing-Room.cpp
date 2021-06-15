#include "LiGu_Graphics\src\RayTracing.h"
#include "LiGu_Graphics\src\Fractal.h"
#pragma GCC optimize(3,"Ofast","inline")
int main() {
	RayTracing ray(1000, 768);
	
	ray.Eye.getData(-500, 0, -50);
	ray.gCenter.getData(0, 0, -50);
	ray.readImg("D:/LIGU.ppm");
	{
		RayTracing::Material* material;
		Mat<double> p1(3, 1), p2(3, 1), p3(3, 1), p4(3, 1), p5(3, 1);

		material = new RayTracing::Material; material->color.getData(1, 1, 1);
		ray.drawSphere(p1.getData(600, -200, -500 + 200), 200, material);
		material = new RayTracing::Material; material->color.getData(1, 1, 1); material->refractRate = 1.6; material->reflectRate = 0.1;
		ray.drawSphere(p1.getData(300, 100, -500 + 100), 100, material);
		material = new RayTracing::Material; material->color.getData(1, 1, 0); material->diffuseReflect = 1;
		ray.drawSphere(p1.getData(500, 300, -500 + 100), 100, material);
		//Box
		material = new RayTracing::Material; material->color.getData(1,0.68,0.75); material->diffuseReflect = 1;
		ray.drawSphere(p1.getData(500, +1e5 + 500, 0), 1E5, material);
		material = new RayTracing::Material; material->color.getData(1, 1, 1); material->diffuseReflect = 1;
		ray.drawSphere(p1.getData(1e5 + 1000, 0, 0), 1E5, material);
		material = new RayTracing::Material; material->color.getData(0.65, 0.62, 1); material->diffuseReflect = 1;
		ray.drawSphere(p1.getData(500, -1e5 - 500, 0), 1E5, material);
		material = new RayTracing::Material; material->color.getData(1, 1, 1); material->diffuseReflect = 1;
		ray.drawSphere(p1.getData(500, 0, +1e5 + 500), 1E5, material);
		ray.drawSphere(p1.getData(500, 0, -1e5 - 500), 1E5, material);
		material = new RayTracing::Material; material->color.getData(10, 10, 10); material->rediateRate = 1;
		ray.drawQuadrilateral(
			p1.getData(500 - 200, -200, 499),
			p2.getData(500 - 200, 200, 499),
			p3.getData(500 + 200, 200, 499),
			p4.getData(500 + 200, -200, 499),
		material);
	}
	ray.paint("D:/LIGU.ppm", 25);
}



		