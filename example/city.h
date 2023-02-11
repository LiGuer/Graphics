#include "C:/UserFiles/Library/Engineering/Graphics/src/Modeling.h"
#include "Academic_Building_Huadian.h"


void city(Modeling& md) {
	Modeling::Point p(3), p1(3);

	double roadL = 3, roadH = 0.25;

	for (int i = 0; i <= 4; i++) {
		for (int j = 0; j <= 4; j++) {
			double
				st_i = 1000.0 * i,
				st_j = 1000.0 * j;

			FourLevelInterchange(md, p = { st_i, st_j, 0 });

			if (j != 4) md.Cuboid(p = { st_i + 2 * roadL, st_j + 500, 0 }, 3 * roadL, 800, roadH);
			if (j != 4) md.Cuboid(p = { st_i - 2 * roadL, st_j + 500, 0 }, 3 * roadL, 800, roadH);
			if (i != 4) md.Cuboid(p = { st_i + 500, st_j + 2 * roadL, 0 }, 800, 3 * roadL, roadH);
			if (i != 4) md.Cuboid(p = { st_i + 500, st_j - 2 * roadL, 0 }, 800, 3 * roadL, roadH);
		}
	}

	double R = 100 - roadL * 2.5;
	vector<Modeling::Point> f;
	f.push_back(p = { R - 2 / 2.0 * roadL,-roadH / 2, 0 });
	f.push_back(p = { R + 2 / 2.0 * roadL,-roadH / 2, 0 });
	f.push_back(p = { R + 2 / 2.0 * roadL, roadH / 2, 0 });
	f.push_back(p = { R - 2 / 2.0 * roadL, roadH / 2, 0 });
	f.push_back(p = { R - 2 / 2.0 * roadL,-roadH / 2, 0 });

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			double
				st_i = 1000.0 * i,
				st_j = 1000.0 * j;

			Grid(md, p = { st_i + 100 ,st_j + 100, 0 }, p1 = { 2 * roadL, 2 * roadL , roadH }, 4, 4, 800, 800);
			md.Rotator(p = { st_i + 500 + R, st_j + 100, 0 }, p1 = { 0, 0, 1 }, f, 720, -PI, -PI / 2);
			md.Rotator(p = { st_i + 500 - R, st_j + 100, 0 }, p1 = { 0, 0, 1 }, f, 720, -PI / 2, 0);
			md.Rotator(p = { st_i + 500 + R, st_j + 900, 0 }, p1 = { 0, 0, 1 }, f, 720, PI / 2, PI);
			md.Rotator(p = { st_i + 500 - R, st_j + 900, 0 }, p1 = { 0, 0, 1 }, f, 720, 0, PI / 2);
			md.Rotator(p = { st_i + 100, st_j + 500 + R, 0 }, p1 = { 0, 0, 1 }, f, 720, PI, PI / 2 * 3);
			md.Rotator(p = { st_i + 100, st_j + 500 - R, 0 }, p1 = { 0, 0, 1 }, f, 720, PI / 2, PI);
			md.Rotator(p = { st_i + 900, st_j + 500 + R, 0 }, p1 = { 0, 0, 1 }, f, 720, -PI / 2, 0);
			md.Rotator(p = { st_i + 900, st_j + 500 - R, 0 }, p1 = { 0, 0, 1 }, f, 720, 0, PI / 2);

			{
				if (i == 1 && j == 1)
					continue;
				vector<double> a(4);

				for (int i_ = 0; i_ < 4; i_++) {
					for (int j_ = 0; j_ < 4; j_++) {
						RandomHouses(md,
							a = {
								st_i + 100 + i_ * 200 + 100,
								st_j + 100 + j_ * 200 + 100,
								190, 190
							},
							20, 40, 10, 3, 15
						);
					}
				}
			}
		}
	}

	// Academic_Building_Huadian_11
	{
		Modeling md_tmp;
		Academic_Building_Huadian_11(md_tmp);

		int i = 1, j = 1, i_ = 3, j_ = 0;
		double
			x = 1000.0 * i + 100 + i_ * 200 + 100 - 120 / 2.0,
			y = 1000.0 * j + 100 + j_ * 200 + 100 - 120 / 2.0;

		md.addTriangleSet(p = { x, y, 0 }, md_tmp.Object);
	}

	// Library_Building_Beiyou_1
	{
		Modeling md_tmp;
		Library_Building_Beiyou_1(md_tmp);

		int i = 1, j = 1;
		double
			x = 1000.0 * i + 500,
			y = 1000.0 * j + 100 + 200 - 3;

		md.addTriangleSet(p = { x, y, 0 }, md_tmp.Object);
	}
}