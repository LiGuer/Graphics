#include "..\LiGu_Codes\LiGu_Graphics\src\RayTracing.h"
#pragma GCC optimize(3,"Ofast","inline")

int main() {
	srand(time(NULL));
	RayTracing ray(700, 700); 
	ray.gCenter = { 100, 100, 100 };
	ray.Eye.add(ray.gCenter, ray.Eye = { -400, 100, 100 });
	{
		Mat<> zero(3), p0(3), p1(3), p2(3), p3(3);
		RayTracing::Material* material;
		material = new RayTracing::Material; material->color = 1; material->reflect = 0.1; material->refractRate = 1.5;
		ray.addCuboid(p1 = { 350, 320,-500 }, p2 = { 650, 350,-200 }, material);
		ray.addCuboid(p1 = { 350,-100,-500 }, p2 = { 650, -70,-200 }, material);
		ray.addCuboid(p1 = { 350,-150,-170 }, p2 = { 650, 400,-200 }, material);
		material = new RayTracing::Material; material->color = { 1, 1, 0 };	material->diffuseReflect = 1;
		ray.addCuboid(p1 = { 400,-400,-500 }, p2 = { 600,-370,-100 }, material);
		ray.addCuboid(p1 = { 400,-200,-500 }, p2 = { 600,-170,-350 }, material);
		ray.addCuboid(p1 = { 400,-400,-350 }, p2 = { 600,-170,-320 }, material);
		material = new RayTracing::Material; material->color = { 1, 1, 1 };	material->reflect = 1;
		ray.addSphere(p1 = { 1000,0,250 }, 200, material);
		//Box
		material = new RayTracing::Material; material->color = { 1, 0.68, 0.75 };	material->diffuseReflect = 1;
		ray.addPlane(p1 = { 0, 1, 0 }, p2 = { 0, +500, 0 }, material);
		material = new RayTracing::Material; material->color = 1;					material->diffuseReflect = 1;
		ray.addPlane(p1 = { 1, 0, 0 }, p2 = { 1000, 0, 0 }, material);
		material = new RayTracing::Material; material->color = { 0.65, 0.62, 1 };	material->diffuseReflect = 1;
		ray.addPlane(p1 = { 0, 1, 0 }, p2 = { 0, -500, 0 }, material);
		material = new RayTracing::Material; material->color = 1;					material->diffuseReflect = 1;
		ray.addPlane(p1 = { 0, 0, 1 }, p2 = { 0, 0, +500 }, material);
		ray.addPlane(p1 = { 0, 0, 1 }, p2 = { 0, 0, -500 }, material);
		material = new RayTracing::Material; material->color = 8; material->rediate = 1;
		ray.addCircle(
			p1 = { 500,0, 499 }, 200,
			p2 = { 0,0,1 },
		material);
		ray.PointLight.push_back(p1 = { 500, 0, 400 });
	}
	ray.paint("D:/LIGU.ppm");
}