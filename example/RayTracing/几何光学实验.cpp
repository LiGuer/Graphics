#include "..\LiGu_Codes\LiGu_Graphics\src\RayTracing.h"
#pragma GCC optimize(3,"Ofast","inline")


// 凸透镜倒立等大的实像
void Model(RayTracing& ray) {
	double a = 600;

	Material* material;
	material = new Material; material->color = 1; material->reflect = 0.1; material->refractRate = 1.5;
	ray.objTree.addStl("D:/tj.stl", { 500, 0, 0 }, 300, &material);
	material = new Material; material->color = { 1, 1, 1 }; material->diffuseReflect = 1;
	ray.objTree.addPlane({ 0, 1, 0 }, { 500, a, 0 }, material);
	material = new Material; material->color = { 1, 0.68, 0.75 }; material->diffuseReflect = 1;
	ray.objTree.addPlane({ 1, 0, 0 }, { 1000, 0, 0 }, material);

	material = new Material; material->color = 20; material->rediate = 1;
	ray.objTree.addCuboid({ 500 - 50, -a - 50, 0 - 50 }, { 500 + 50, -a + 50, 0 + 50 }, material);

	material = new Material; material->color = { 0, 0, 20 }; material->rediate = 1;
	ray.objTree.addSphere({ 500, -a, 150 }, 50, material);

	material = new Material; material->color = { 0, 20, 0 }; material->rediate = 1;
	ray.objTree.addSphere({ 500 - 100, -a, -100 }, 50, material);

	material = new Material; material->color = { 20, 0, 0 }; material->rediate = 1;
	ray.objTree.addSphere({ 500 + 100, -a, -100 }, 50, material);

}

int main() {
	srand(time(NULL));
	RayTracing ray(800, 500);
	ray.gCenter = { 0, 0, 0 };
	ray.Eye.add(ray.gCenter, ray.Eye = { -400, 0, 0 });
	Model(ray);

	ray.paint("D:/LIGU_ROOM.ppm");
}