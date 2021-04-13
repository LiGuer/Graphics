#include"src/RayTracing.h"

int main() {
	RayTracing ray(1000, 1000);
	
	{ double t[] = { -600,0,0 }; ray.Eye.getData(t); }
	{ double t[] = { -100,0,0 }; ray.gCenter.getData(t); }

	{
		RayTracing::Material* material;
		Mat<double> p1(3, 1), p2(3, 1), p3(3, 1), p4(3, 1), p5(3, 1);

		material = new RayTracing::Material;
		material->color = 0xFFFFCC; 
		{ double t[] = { 600,-200,-500+200 }; p1.getData(t); }
		ray.drawSphere(p1, 200, material);

		material = new RayTracing::Material;
		material->color = 0xFFFFFF; material->refractRate = 1.6; material->reflectRate = 0.1;
		{ double t[] = { 200,100,-500+100 }; p1.getData(t); }
		ray.drawSphere(p1, 100, material);


		material = new RayTracing::Material;
		material->color = 0x00FF00;
		{ double t[] = { 300,-300,-500 }; p1.getData(t); }
		{ double t[] = { 400,-200,-400 }; p2.getData(t); }
		//ray.drawCuboid(p1, p2, material);

		material = new RayTracing::Material;
		material->color = 0xFFCCCC; material->refractRate = 1.6;  material->reflectRate = 0.3;
		{ double t[] = { 500,-500,-500 }; p1.getData(t); }
		{ double t[] = { 600,-400,-400 }; p2.getData(t); }
		//ray.drawCuboid(p1, p2, material);


		material = new RayTracing::Material;
		material->color = 0xFFCCFF; material->refractRate = 1.6; material->reflectRate = 0;
		{ double t[] = { 800,500,-500 + 100 }; p1.getData(t); }
		//ray.drawSphere(p1, 100, material);


		material = new RayTracing::Material;
		material->color = 0x8cfffb; 
		{ double t[] = { 900,800,-500 + 100 }; p1.getData(t); }
		//ray.drawSphere(p1, 100, material);
		

		material = new RayTracing::Material;
		material->color = 0xFFFFEE; 
		{ double t[] = { 0,-500,-500 }; p1.getData(t); }
		{ double t[] = { 0,500,-500 }; p2.getData(t); }
		{ double t[] = { 1000,500,-500 }; p3.getData(t); }
		{ double t[] = { 1000,-500,-500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);


		material = new RayTracing::Material;
		material->color = 0xFFFFFF; material->diffuseReflect = 1;
		{ double t[] = { 0,-500,500 }; p1.getData(t); }
		{ double t[] = { 0,500,500 }; p2.getData(t); }
		{ double t[] = { 1000,500,500 }; p3.getData(t); }
		{ double t[] = { 1000,-500,500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);

		material = new RayTracing::Material;
		material->color = 0xFF0000; material->diffuseReflect = 1;
		{ double t[] = { 0,-500,-500 }; p1.getData(t); }
		{ double t[] = { 0,-500,500 }; p2.getData(t); }
		{ double t[] = { 1000,-500,500 }; p3.getData(t); }
		{ double t[] = { 1000,-500,-500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);

		material = new RayTracing::Material;
		material->color = 0x0000FF; material->diffuseReflect = 1;
		{ double t[] = { 0,500,-500 }; p1.getData(t); }
		{ double t[] = { 0,500,500 }; p2.getData(t); }
		{ double t[] = { 1000,500,500 }; p3.getData(t); }
		{ double t[] = { 1000,500,-500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);

		material = new RayTracing::Material;
		material->color = 0x00EE00; material->diffuseReflect = 1;
		{ double t[] = { 1000,-500,-500 }; p1.getData(t); }
		{ double t[] = { 1000,500,-500 }; p2.getData(t); }
		{ double t[] = { 1000,500,500 }; p3.getData(t); }
		{ double t[] = { 1000,-500,500 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);



		material = new RayTracing::Material;
		material->color = 0xFFFFFF; material->rediateRate = 1;
		{ double t[] = { 500 - 100,-100,499 }; p1.getData(t); }
		{ double t[] = { 500 - 100,100,499 }; p2.getData(t); }
		{ double t[] = { 500 + 100,100,499 }; p3.getData(t); }
		{ double t[] = { 500 + 100,-100,499 }; p4.getData(t); }
		ray.drawQuadrilateral(p1, p2, p3, p4, material);



		{ double t[] = { 500,0,500 }; p1.getData(t); }
		ray.LightSource.push_back(p1);

	}

	ray.paint();
	ray.g.writeImg("D:/LIGU.ppm");

}