#ifndef STAIR_H
#define STAIR_H

#include "../src/Modeling.h"

inline void Stair(Modeling& md,
	Modeling::Point& certer, vector<double>& direction,
	double L, double W, double H, int numStep
) {
	Modeling::Point p(3), delta(3);
	double
		stapL = L / numStep,
		stapH = H / numStep,
		theta = (direction[1] > 0 ? 1 : -1) * acos(direction[0] / sqrt(pow(direction[0], 2) + pow(direction[1], 2)));

	delta = {
		cos(theta) * stapL,
		sin(theta) * stapL,
		stapH
	};
	add(p, certer, mul(p, 0.5, delta));

	for (int i = 0; i < numStep; i++) {
		md.Cuboid(p, direction, stapL, W, stapH);
		add(p, p, delta);
	}
}

inline void Stair_2(Modeling& md, Modeling::Point certer, vector<double> direction,
	double L, double W, double H, int numLayer, int numStep,
	double platformL, int leftfirst = 1
) {
	if (leftfirst != 1)
		leftfirst = -1;

	Modeling::Point p(3), p1(3), p2(3), p3(3);
	Mat<double> rotateMat(3, 3);
	double
		L_ = L - platformL,
		theta = (direction[1] > 0 ? 1 : -1) * acos(direction[0] / sqrt(pow(direction[0], 2) + pow(direction[1], 2)));

	rotateMat = {
		cos(theta),-sin(theta), 0,
		sin(theta), cos(theta), 0,
		0, 0, 1
	};

	for (int k = 0; k < numLayer; k++) {
		mul(p, rotateMat, p = { 0, leftfirst * W / 4, 0 });
		add(p, certer, p);
		Stair(md, p, direction, L_, W / 2, H / 2, numStep / 2);

		mul(p, rotateMat, p = { L_ + platformL / 2 ,0, H / 2 });
		add(p, certer, p);
		md.Cuboid(p, direction, platformL, W, H / numStep);

		mul(direction, -1, direction);

		mul(p, rotateMat, p = { L_, -leftfirst * W / 4, H / 2 });
		add(p, certer, p);
		Stair(md, p, direction, L_, W / 2, H / 2, numStep / 2);

		mul(direction, -1, direction);
		certer[2] += H;
	}

}

inline void Stair_3(Modeling& md, Modeling::Point certer, vector<double> direction,
	double L, double W, double H, int numLayer, int numStep,
	double platformL
) {
	Modeling::Point p(3), p1(3), p2(3), p3(3);
	Mat<double> rotateMat(3, 3);
	double
		L_ = L - platformL,
		theta = (direction[1] > 0 ? 1 : -1) * acos(direction[0] / sqrt(pow(direction[0], 2) + pow(direction[1], 2)));

	rotateMat = {
		cos(theta),-sin(theta), 0,
		sin(theta), cos(theta), 0,
		0, 0, 1
	};

	for (int k = 0; k < numLayer; k++) {
		Stair(md, p = certer, direction, L_, W / 2, H / 2, numStep / 2);

		mul(p, rotateMat, p = { L_ + platformL / 2 ,0, H / 2 });
		add(p, certer, p);
		md.Cuboid(p, direction, platformL, W, H / numStep);

		mul(direction, -1, direction);

		mul(p, rotateMat, p = { L_, W / 4 * 1.5, H / 2 });
		add(p, certer, p);
		Stair(md, p, direction, L_, W / 4, H / 2, numStep / 2);


		mul(p, rotateMat, p = { L_, -W / 4 * 1.5, H / 2 });
		add(p, certer, p);
		Stair(md, p, direction, L_, W / 4, H / 2, numStep / 2);

		mul(direction, -1, direction);
		certer[2] += H;
	}
}

inline void Stair_Spiral(Modeling& md, Modeling::Point certer,
	double R_in, double R_out, double H,
	double angleSt, double angleEd, int numStep
) {
	Modeling::Point p(3), p1(3), p2(3), p3(3);
	vector<Modeling::Point> f;
	Mat<double> rotateMat(3, 3);

	double h = H / numStep;

	f.push_back({ R_in, -h / 2, 0 });
	f.push_back({ R_out, -h / 2, 0 });
	f.push_back({ R_out,  h / 2, 0 });
	f.push_back({ R_in,  h / 2, 0 });
	f.push_back({ R_in, -h / 2, 0 });

	for (int i = 0; i < numStep; i++) {
		certer[2] += h;

		double
			x = (double)i / numStep,
			a = x * (angleEd - angleSt) + angleSt;

		md.Rotator(certer, p = { 0, 0, 1 }, f, 60,
			a, a + (angleEd - angleSt) / numStep, true);
	}
}

#endif // !STAIR_H
