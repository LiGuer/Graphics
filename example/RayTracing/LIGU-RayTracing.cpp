#include"src/RayTracing.h"

int main() {
	RayTracing ray(1000, 1000);
	
	{ double t[] = { -500,0,0 }; ray.Eye.getData(t); }
	{ double t[] = { 0,0,0 }; ray.gCenter.getData(t); }

	{
		RayTracing::Material* material;
		Mat<double> p1(3, 1), p2(3, 1), p3(3, 1), p4(3, 1), p5(3, 1);

		material = new RayTracing::Material;
		material->color = 0xFFFFCC; 
		{ double t[] = { 500,0,-500+200 }; p1.getData(t); }
		ray.drawSphere(p1, 200, material);

		material = new RayTracing::Material;
		material->color = 0xFFFFFF; material->refractiveIndex = 1.6; material->reflectance = 0.1;
		{ double t[] = { 300,500,-500+100 }; p1.getData(t); }
		ray.drawSphere(p1, 100, material);


		material = new RayTracing::Material;
		material->color = 0x00FF00;
		{ double t[] = { 300,-300,-500 }; p1.getData(t); }
		{ double t[] = { 400,-200,-400 }; p2.getData(t); }
		ray.drawCuboid(p1, p2, material);

		material = new RayTracing::Material;
		material->color = 0xFFCCCC; material->refractiveIndex = 1.6;  material->reflectance = 0.3;
		{ double t[] = { 500,-500,-500 }; p1.getData(t); }
		{ double t[] = { 600,-400,-400 }; p2.getData(t); }
		ray.drawCuboid(p1, p2, material);


		material = new RayTracing::Material;
		material->color = 0xFFCCFF; material->refractiveIndex = 1.6; material->reflectance = 0;
		{ double t[] = { 800,500,-500 + 100 }; p1.getData(t); }
		ray.drawSphere(p1, 100, material);


		material = new RayTracing::Material;
		material->color = 0x8cfffb; 
		{ double t[] = { 900,800,-500 + 100 }; p1.getData(t); }
		ray.drawSphere(p1, 100, material);
		

		material = new RayTracing::Material;
		material->color = 0xFFFFEE; 
		{ double t[] = { -3000,-5000,-500 }; p1.getData(t); }
		{ double t[] = { -3000,5000,-500 }; p2.getData(t); }
		{ double t[] = { 5000,5000,-500 }; p3.getData(t); }
		{ double t[] = { 5000,-5000,-500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);



		{ double t[] = { -800,-1000,5000 }; p1.getData(t); }
		ray.LightSource.push_back(p1);

	}

	ray.paint();
	ray.g.writeImg("D:/LIGU.ppm");

}