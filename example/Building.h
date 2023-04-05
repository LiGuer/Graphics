#include "C:/UserFiles/Library/Engineering/Graphics/src/Modeling.h"
#include "Stair.h"


void pipe(Modeling& md, vector<double>& st, vector<double>& ed, double r_in, double r_out) {
	vector<vector<double>> f, path;

	md.Circle(f, r_out, 720);
	f.push_back({ 0, r_in });
	md.Circle(f, r_in, 720, 1);
	f.push_back({ 0, r_out });

	path.push_back(st);
	path.push_back(ed);

	md.Translator(path, f);
}


inline void FourLevelInterchange(Modeling& md, Modeling::Point& c, double roadL = 3, double roadH = 0.25, double tatalL = 200, double layerH = 8) {
	Modeling::Point p1(3), p2(3), p3(3);
	double
		x_step[] = { 1, 1,-1,-1 },
		y_step[] = { 1,-1, 1,-1 };
	int N = 720;

	std::function<double(double)> func = [](double x) {
		return exp(-x) / pow(1 + exp(-x), 2);
	};

	// Straight
	{
		vector<Modeling::Point> f;

		f.push_back(p1 = { -3 / 2.0 * roadL,-roadH / 2, 0 });
		f.push_back(p1 = { +3 / 2.0 * roadL,-roadH / 2, 0 });
		f.push_back(p1 = { +3 / 2.0 * roadL, roadH / 2, 0 });
		f.push_back(p1 = { -3 / 2.0 * roadL, roadH / 2, 0 });
		f.push_back(p1 = { -3 / 2.0 * roadL,-roadH / 2, 0 });

		int N = 720;

		for (int k = 0; k < 4; k++) {
			vector<Modeling::Point> path;

			for (int i = 0; i < N; i++) {
				double
					a = (i / (double)N - 0.5) * 2,
					x = a * 8;

				p1 = {
					c[0],
					c[1] + a * tatalL / 2.0,
					c[2]
				};

				if (k == 1 || k == 3)
					p1[0] += 0.5 * roadL + 3 / 2.0 * roadL;
				else
					p1[0] -= 0.5 * roadL + 3 / 2.0 * roadL;

				if (k >= 2) {
					swap(p1[0], p1[1]);
					p1[2] += 3 * layerH * func(x);
				}

				path.push_back(p1);
			}
			md.Translator(path, f);
		}
	}

	// Right
	{
		vector<Modeling::Point> f;

		f.push_back(p1 = { -2 / 2.0 * roadL,-roadH / 2, 0 });
		f.push_back(p1 = { +2 / 2.0 * roadL,-roadH / 2, 0 });
		f.push_back(p1 = { +2 / 2.0 * roadL, roadH / 2, 0 });
		f.push_back(p1 = { -2 / 2.0 * roadL, roadH / 2, 0 });
		f.push_back(p1 = { -2 / 2.0 * roadL,-roadH / 2, 0 });

		double R = tatalL / 2.0;

		for (int k = 0; k < 4; k++) {
			vector<Modeling::Point> path;

			for (int i = 0; i < N; i++) {
				double a = 0.025 * PI + PI * 0.45 / (double)N * i;

				p1 = {
					c[0] + x_step[k] * (1.5 * roadL + R - R * sin(a)),
					c[1] + y_step[k] * (1.5 * roadL + R - R * cos(a)),
					c[2]
				};

				path.push_back(p1);
			}
			md.Translator(path, f);
		}
	}

	// Left
	{
		vector<Modeling::Point> f;
		f.push_back(p1 = { -1 / 2.0 * roadL,-roadH / 2, 0 });
		f.push_back(p1 = { +1 / 2.0 * roadL,-roadH / 2, 0 });
		f.push_back(p1 = { +1 / 2.0 * roadL, roadH / 2, 0 });
		f.push_back(p1 = { -1 / 2.0 * roadL, roadH / 2, 0 });
		f.push_back(p1 = { -1 / 2.0 * roadL,-roadH / 2, 0 });

		double R = tatalL * 0.3;

		for (int k = 0; k < 4; k++) {
			vector<Modeling::Point> path;

			for (int i = 0; i < N; i++) {
				double
					a = PI / 4 + (i - N / 2) / (double)N * 2.2,
					x = (i / (double)N - 0.5) * 2 * 6;

				p1 = {
					c[0] + x_step[k] * (roadL / sqrt(2) + R / sqrt(2) - R * sin(a)),
					c[1] + y_step[k] * (roadL / sqrt(2) + R / sqrt(2) - R * cos(a)),
					c[2]
				};

				if (k == 1 || k == 2)
					p1[2] += 2 * layerH * func(x);
				else
					p1[2] += layerH * func(x);

				path.push_back(p1);
			}
			md.Translator(path, f);
		}
	}
}


void Grid(Modeling& md, Modeling::Point& st, vector<double>& edgeLen, int nx, int ny, double len_x, double len_y, vector<double>& direction) {
	Modeling::Point p(3);
	double
		gridLen_x = len_x / nx,
		gridLen_y = len_y / ny;

	for (int i = 0; i <= nx; i++)
		md.Cuboid(
			p = { st[0] + i * gridLen_x, st[1] + len_y / 2, st[2] },
			edgeLen[0], len_y + edgeLen[1], edgeLen[2]
		);

	for (int i = 0; i <= ny; i++)
		md.Cuboid(
			p = { st[0] + len_y / 2, st[1] + i * gridLen_y, st[2] }, 
			len_x + edgeLen[0], edgeLen[1], edgeLen[2]
		);
}

// Random Houses
void RandomHouses(
	Modeling& md, vector<double>& rectangle, 
	double maxHeight, double minHeight, int N,
	double roadsize = 0, int sizeThreshold = 0, int randThreshold = 0
) {
	vector<double> p(4), p1(4), p2(4), t1(3), t2(3);
	vector<vector<double>> rectangleSet;
	rectangleSet.push_back(rectangle);

	for (int i = 0; i < N; i++) {
		int n = rectangleSet.size();

		for (int j = 0; j < n; j++) {
			p1 = p2 = rectangleSet[j];

			int ind = (i % 2 == 0) ? 0 : 1;

			if (rectangleSet[j][ind + 2] - roadsize > sizeThreshold) {
				double ra = rand() / (double)RAND_MAX * (1 - 2 * randThreshold) + randThreshold;
				ra *= rectangleSet[j][ind + 2];

				p1[ind + 2] = ra;
				p2[ind + 2] = rectangleSet[j][ind + 2] - ra;
				p1[ind] = rectangleSet[j][ind] - rectangleSet[j][ind + 2] / 2 + ra / 2;
				p2[ind] = rectangleSet[j][ind] + ra / 2;

				rectangleSet[j] = p1;
				rectangleSet.push_back(p2);
			}
		}
	}

	int n = rectangleSet.size();
	for (int j = 0; j < n; j++) {
		if (rectangleSet[j][2] - roadsize < sizeThreshold / 2 ||
			rectangleSet[j][3] - roadsize < sizeThreshold / 2)
			continue;

		double h = rand() / (double)RAND_MAX * (maxHeight - minHeight) + minHeight;

		md.Cuboid(
			t1 = { rectangleSet[j][0], rectangleSet[j][1], h / 2 },
			rectangleSet[j][2] - roadsize, rectangleSet[j][3] - roadsize, h
		);
	}
}