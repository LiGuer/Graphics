#include "..\LiGu_Codes\LiGu_Graphics\src\RayTracing.h"
#pragma GCC optimize(3,"Ofast","inline")
int main() {
	srand(time(NULL));
	RayTracing ray(1,1);
	ray.Eye = { -400, 0, -50 };
	ray.gCenter = { 0, 0, -50 };
	GraphicsFileCode::ppmRead("D:/LIGU.ppm", ray.ScreenPix); ray.Screen.alloc(ray.ScreenPix.rows, ray.ScreenPix.cols);
	for (int i = 0; i < ray.ScreenPix.size(); i++)
		ray.Screen(ray.Screen.rows - 1 - ray.Screen.i2x(i), ray.Screen.i2y(i)).zero(3) = { ray.ScreenPix[i].R / 255.0, ray.ScreenPix[i].G / 255.0, ray.ScreenPix[i].B / 255.0 };
	{
		RayTracing::Material* material;
		Mat<> p1(3), p2(3), p3(3), p4(3), p5(3);

		material = new RayTracing::Material; material->color = 1;
		ray.drawSphere(p1 = { 600, -200, -500 + 200 }, 200, material);
		material = new RayTracing::Material; material->color = 1;					material->refractRate = 1.6; material->reflect = 0.1;
		ray.drawSphere(p1 = { 300, 100, -500 + 100 }, 100, material);
		material = new RayTracing::Material; material->color = { 1, 1, 0 };			material->diffuseReflect = 1;
		ray.drawSphere(p1 = { 500, 300, -500 + 100 }, 100, material);
		//Box
		material = new RayTracing::Material; material->color = { 1, 0.68, 0.75 };	material->diffuseReflect = 1;
		ray.drawSphere(p1 = { 500, +1e5 + 500, 0 }, 1E5, material);
		material = new RayTracing::Material; material->color = 1;					material->diffuseReflect = 1;
		ray.drawSphere(p1 = { 1e5 + 1000, 0, 0 }, 1E5, material);
		material = new RayTracing::Material; material->color = { 0.65, 0.62, 1 };	material->diffuseReflect = 1;
		ray.drawSphere(p1 = { 500, -1e5 - 500, 0 }, 1E5, material);
		material = new RayTracing::Material; material->color = 1;					material->diffuseReflect = 1;
		ray.drawSphere(p1 = { 500, 0, +1e5 + 500 }, 1E5, material);
		ray.drawSphere(p1 = { 500, 0, -1e5 - 500 }, 1E5, material);
		material = new RayTracing::Material; material->color = 12;					material->rediate = 1;
		ray.drawTriangle(
			p1 = { 500 - 200, -200, 499 },
			p2 = { 500 - 200, 200, 499 },
			p3 = { 500 + 200, 200, 499 },
			material);
	}
	ray.paint("D:/LIGU.ppm", 2630);
}