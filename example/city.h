#include "C:/UserFiles/Library/Engineering/Graphics/src/Modeling.h"
#include "Academic_Building_Huadian.h"


void city(Modeling& md) {
	Modeling::Point p(3), p1(3), p2(3);
	int _N = 1;

	{// freeway
		double H = 16, roadL = 3, roadH = 0.25;

		for (int i = 0; i <= _N; i++) {
			for (int j = 0; j <= _N; j++) {
				double
					st_i = 1000.0 * i,
					st_j = 1000.0 * j;

				FourLevelInterchange(md, p = { st_i, st_j, H });

				if (j != _N) md.Cuboid(p = { st_i + 2 * roadL, st_j + 500, H }, 3 * roadL, 800, roadH);
				if (j != _N) md.Cuboid(p = { st_i - 2 * roadL, st_j + 500, H }, 3 * roadL, 800, roadH);
				if (i != _N) md.Cuboid(p = { st_i + 500, st_j + 2 * roadL, H }, 800, 3 * roadL, roadH);
				if (i != _N) md.Cuboid(p = { st_i + 500, st_j - 2 * roadL, H }, 800, 3 * roadL, roadH);
			}
		}

		// gateway
	}

	{//center roads
		double roadL = 3, roadH = 0.25;

		for (int i = 0; i < _N; i++)
			md.Cuboid(p = { 1000.0 * i + 500, 1000 * _N / 2.0, roadH / 2 }, 2 * roadL, 1000 * _N, roadH);
		for (int i = 0; i < _N; i++)
			md.Cuboid(p = { 1000 * _N / 2.0, 1000.0 * i + 500, roadH / 2 }, 1000 * _N, 2 * roadL, roadH);
	}

	{
		// subway
		double
			H = -32, dH = -8, h = 6,
			R = h / 2, dR = 0.3;
		vector<Modeling::Point> f_out, f_in;
		{
			Modeling::Circle(f_in,  R, 120);
			Modeling::Circle(f_out, R - dR, 120);
		}

		for (int i = 0; i < _N; i++) {
			for (int j = 0; j < _N; j++) {
				double
					st_i = 1000.0 * i + 500,
					st_j = 1000.0 * j + 500,
					step_i[] = {500,   0, 500, -500},
					step_j[] = {  0, 500, 500,  500},
					delta = R + 12 / 2.0,
					step_i_2[] = { 0, delta,-delta / sqrt(2), delta / sqrt(2) },
					step_j_2[] = { delta, 0, delta / sqrt(2), delta / sqrt(2) };
					

				for (int k = 0; k < 8; k++) {
					double
						st_x = st_i - step_i[k % 4] + step_i_2[k % 4] * (k < 4 ? 1 : -1),
						ed_x = st_i + step_i[k % 4] + step_i_2[k % 4] * (k < 4 ? 1 : -1),
						st_y = st_j - step_j[k % 4] + step_j_2[k % 4] * (k < 4 ? 1 : -1),
						ed_y = st_j + step_j[k % 4] + step_j_2[k % 4] * (k < 4 ? 1 : -1),
						h_t = H + (k % 4) * dH;

					md.Translator(
						p1 = { st_x, st_y, h_t },
						p2 = { ed_x, ed_y, h_t }, f_in, 0);
					md.Translator(
						p1 = { st_x, st_y, h_t },
						p2 = { ed_x, ed_y, h_t }, f_out, 0);
				}
			}
		}

		// subway station
		{
			double h = 0.25,
				step_i[] = { 1, 0, 1,-1 },
				step_j[] = { 0, 1, 1, 1 };

			for (int i = 0; i < _N; i++) {
				for (int j = 0; j < _N; j++) {
					double
						st_i = 1000.0 * i + 500,
						st_j = 1000.0 * j + 500;

					md.Cuboid(p1 = { st_i, st_j, H - dH - h / 2 - dH / 2.0 }, 64, 64, 6);
					md.Cuboid(p1 = { st_i, st_j, H - dH - h / 2 - dH / 2.0 }, 64, 64, 6 - 2 * h);

					for (int k = 0; k < 4; k++) {
						md.Cuboid(p1 = { st_i, st_j, H - R + k * dH }, p = { step_i[k], step_j[k], 0 }, 256, 12 + 2 * R, h);
						md.Cuboid(p1 = { st_i, st_j, H + R + k * dH }, p = { step_i[k], step_j[k], 0 }, 256, 12 + 2 * R, h);
					}
					
				}
			}
		}
	}

	{// Random Houses
		for (int i = 0; i < _N; i++) {
			for (int j = 0; j < _N; j++) {
				double
					st_i = 1000.0 * i,
					st_j = 1000.0 * j;
				vector<double> a(4);
				/*
				RandomHouses(md,
					a = {
						st_i + 500,
						st_j + 500,
						800, 800
					},
					20, 80, 20, 3, 15
				);*/

			}
		}
	}
}