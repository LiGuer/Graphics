#include "..\LiGu_Codes\LiGu_Graphics\src\RayTracing.h"
#pragma GCC optimize(3,"Ofast","inline")
Mat<bool> img1, img2, img3, img4;

bool f1(double x, double y) { std::swap(x, y); x *= 1.5; y *= 1.5; if (x < 0 || y < 0 || x >= img1.rows || y >= img1.cols) return 0; return img1(img1.rows - 1 - x, y); }
bool f3(double x, double y) { if (x < PI / 6)return true; if (x > PI - PI / 6) return 0; if (y > PI) y -= PI; return img3(y * (img3.rows - 1) / PI, (img3.cols - 1) / (2 * PI / 3) * (x - PI / 6)); }

void Model(RayTracing& ray) {
	Material* material;
	material = new Material; material->color = 1;					material->reflect = 0.1; material->refractRate = 1.5;
	ray.objTree.addCuboid({ 350, 320,-500 }, { 650, 350,-200 }, material);
	ray.objTree.addCuboid({ 350,-100,-500 }, { 650, -70,-200 }, material);
	ray.objTree.addCuboid({ 350,-150,-200 }, { 650, 400,-170 }, material);
	ray.objTree.addPlane({ 1, 0, 0 }, { 1000, 0, 0 }, material);
	ray.objTree.addPlane({ 1, 0, 0 }, { 950, 0, 0 }, material);
	ray.objTree.addStl("D:/bun_zipper2.stl", { 300, 150,-500 }, 3, &material);
	material = new Material; material->color = { 1, 1, 0 };			material->diffuseReflect = 1;
	ray.objTree.addCuboid({ 400,-400,-500 }, { 600,-370,-100 }, material);
	ray.objTree.addCuboid({ 400,-200,-500 }, { 600,-170,-350 }, material);
	ray.objTree.addCuboid({ 400,-400,-350 }, { 600,-170,-320 }, material);
	material = new Material; material->color = { 1, 1, 0 };			material->reflect = 1;
	ray.objTree.addPlaneShape({ 0, 1, 0 }, { 100, -495, -250 }, f1, material);
	material = new Material; material->color = { 1, 1, 1 };			material->reflect = 1;
	ray.objTree.addSphere({ 500, 500 + 75, 200 }, 150, material);
	material = new Material; material->color = { 1, 0.678, 0.784 };	material->reflect = 1;
	ray.objTree.addSphere({ 500,250,-170 + 120 }, 120, material, [](double x, double y) {
		double o = y - ((int)(y / (PI / 2))) * (PI / 2), r = 3 / 2.5 * sin(o) * cos(o) / (pow(sin(o), 3) + pow(cos(o), 3)); if ((PI - x) < asin(r)) return true; return false; });
	ray.objTree.addSphere({ 500,250,-170 + 120 / 3.0 }, 120 / 3.0, material, [](double x, double y) { y += PI / 4;
	double o = y - ((int)(y / (PI / 2))) * (PI / 2), r = 3 / 2.5 * sin(o) * cos(o) / (pow(sin(o), 3) + pow(cos(o), 3)); if ((PI - x) < 2 * asin(r)) return true; return false; });
	material = new Material; material->color = 1;	material->rediate = 1;
	ray.objTree.addSphere({ 500,250,-20 }, 30, material);
	material = new Material; material->color = { 1, 0.678, 0.784 };	material->diffuseReflect = 1;
	ray.objTree.addSphere({ 500,250,-20 }, 70, material, f3);
	{
		Mat<RGB> img0; GraphicsFileCode::ppmRead("D:/f_img.ppm", img0);
		img1.zero(img0.rows, img0.cols); for (int i = 0; i < img1.size(); i++) img1[i] = img0[i].R < 100 ? 1 : 0;

		GraphicsFileCode::ppmRead("D:/f_img3.ppm", img0); img0.transpose(img0);
		for (int x = 0; x < img0.rows; x++) for (int y = 0; y < img0.cols / 2; y++) std::swap(img0(x, y), img0(x, img0.cols - 1 - y));
		img3.zero(img0.rows, img0.cols); for (int i = 0; i < img3.size(); i++) img3[i] = img0[i].R < 100 ? 1 : 0;
	}
	material = new Material; material->color = { 1, 1, 1 };	material->diffuseReflect = 1;
	ray.objTree.addStl("D:/teapot.stl", { 500, 0,-170 }, 10, &material);
	ray.objTree.addStl("D:/Venus2.stl", { 600, -400,-500 }, 3.2, &material);
	//Box
	material = new Material; material->color = { 0.65, 0.62, 1 };	material->diffuseReflect = 1;
	ray.objTree.addPlaneShape({ 0, 1, 0 }, { 0, +500, 0 }, [](double x, double y) { return x <= 1000 ? true : false; }, material);
	material = new Material; material->color = { 1, 0.68, 0.75 };	material->diffuseReflect = 1;
	ray.objTree.addPlaneShape({ 0, 1, 0 }, { 0, -500, 0 }, [](double x, double y) { return x <= 1000 ? true : false; }, material);
	material = new Material; material->color = 1;					material->diffuseReflect = 1;
	ray.objTree.addPlane({ 1, 0, 0 }, { 2000, 0, 0 }, material);
	ray.objTree.addPlane({ 1, 0, 0 }, { -100, 0, 0 }, material);
	ray.objTree.addPlane({ 0, 0, 1 }, { 0, 0, +500 }, material);
	ray.objTree.addPlane({ 0, 0, 1 }, { 0, 0, -500 }, material);
	material = new Material; material->color = 7;					material->rediate = 1;
	ray.objTree.addCircle({ 500,0, 499 }, 250, { 0,0,1 }, material);
	ray.objTree.addCircle({ 1500,0, 499 }, 250, { 0,0,1 }, material);

	//Mat<> p(3);
	//ray.PointLight.push_back(p = { 500, 0, 400 });
}

int main() {
	srand(time(NULL));
	RayTracing ray(800, 800);
	ray.gCenter = { 100, 100, 100 };
	ray.Eye.add(ray.gCenter, ray.Eye = { -400, 100, 100 });
	Model(ray);
	ray.paint("D:/LIGU_ROOM.ppm");
}