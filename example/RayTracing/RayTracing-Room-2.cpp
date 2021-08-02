#include "..\LiGu_Codes\LiGu_Graphics\src\RayTracing.h"
#pragma GCC optimize(3,"Ofast","inline")
Mat<bool> img1, img2;
bool f1(double x, double y) {
	if (x * 1.7 < 0 || y * 1.7 < 0) return 0;
	if (x * 1.7 >= img1.rows || y * 1.7 >= img1.cols) return 0;
	return img1(x * 1.7, y * 1.7);
}
bool f2(double x, double y) {
	if (x < 0 || y < 0) return 0;
	if (x >= img2.rows || y >= img2.cols) return 0;
	return img2(x, y);
}

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
		ray.addSphere(p1 = { 500, 500 + 75,200 }, 150, material);
		material = new RayTracing::Material; material->color = { 1, 0.678, 0.784 };	material->reflect = 1;
		ray.addSphere(p1 = { 500,150,-200 + 300 }, 300, material, [](double x, double y) {
			double o = y - ((int)(y / (PI / 2))) * (PI / 2);
			double r = 3 / 2.5 * sin(o) * cos(o) / (pow(sin(o), 3) + pow(cos(o), 3));
			if ((PI - x) < asin(r)) return true; return false;
			});
		ray.addSphere(p1 = { 500,150,-200 + 100 }, 100, material, [](double x, double y) {
			y += PI / 4;
			double o = y - ((int)(y / (PI / 2))) * (PI / 2);
			double r = 3 / 2.5 * sin(o) * cos(o) / (pow(sin(o), 3) + pow(cos(o), 3));
			if ((PI - x) < 2 * asin(r)) return true; return false;
			});
		{
			Mat<RGB> img0;
			GraphicsFileCode::ppmRead("D:/f_img.ppm", img0);
			img1.zero(img0.rows, img0.cols); for (int i = 0; i < img1.size(); i++) img1[i] = img0[i].R < 100 ? 1 : 0;
			material = new RayTracing::Material; material->color = { 1, 1, 0 };	material->reflect = 1;
			ray.addPlaneShape(p1 = { 0, -495, -500 }, p2 = { 0, 1, 0 }, f1, material);

			GraphicsFileCode::ppmRead("D:/f_img2.ppm", img0); img0.transpose(img0);
			for (int x = 0; x < img0.rows; x++)
				for (int y = 0; y < img0.cols / 2; y++)
					std::swap(img0(x, y), img0(x, img0.cols - 1 - y));
			img2.zero(img0.rows, img0.cols); for (int i = 0; i < img2.size(); i++) img2[i] = img0[i].R < 100 ? 1 : 0;
			material = new RayTracing::Material; material->color = { 1, 0.678, 0.784 };	material->reflect = 1;
			ray.addPlaneShape(p1 = { 995, -450, -200 }, p2 = { -1, 0, 0 }, f2, material);
		}
		//Box
		material = new RayTracing::Material; material->color = { 1, 0.68, 0.75 };	material->diffuseReflect = 1;
		ray.addPlane(p1 = { 0, 1, 0 }, p2 = { 0, -500, 0 }, material);
		material = new RayTracing::Material; material->color = 1;					material->diffuseReflect = 1;
		ray.addPlane(p1 = { 1, 0, 0 }, p2 = { 1000, 0, 0 }, material);
		ray.addPlane(p1 = { 1, 0, 0 }, p2 = { -100, 0, 0 }, material);
		material = new RayTracing::Material; material->color = { 0.65, 0.62, 1 };	material->diffuseReflect = 1;
		ray.addPlane(p1 = { 0, 1, 0 }, p2 = { 0, +500, 0 }, material);
		material = new RayTracing::Material; material->color = 1;					material->diffuseReflect = 1;
		ray.addPlane(p1 = { 0, 0, 1 }, p2 = { 0, 0, +500 }, material);
		ray.addPlane(p1 = { 0, 0, 1 }, p2 = { 0, 0, -500 }, material);
		material = new RayTracing::Material; material->color = 7; material->rediate = 1;
		ray.addCircle(
			p1 = { 500,0, 499 }, 250,
			p2 = { 0,0,1 },
			material);
		ray.PointLight.push_back(p1 = { 500, 0, 400 });
	}
	ray.paint("D:/LIGU.ppm");
}