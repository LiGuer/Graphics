#include "C:/UserFiles/Library/Engineering/Graphics/src/Modeling.h"
#include "Building.h"

void Academic_Building_Huadian_11(Modeling& md) {
	vector<Modeling::Point> f;
	Modeling::Point p(3), p1(3), p2(3);

	// wall
	{
		double wallW = 0.2;

		f.push_back(p = { -1 / 2.0 * wallW,0, 0 });
		f.push_back(p = { 1 / 2.0 * wallW,0, 0 });
		f.push_back(p = { 1 / 2.0 * wallW, 1.2, 0 });
		f.push_back(p = { -1 / 2.0 * wallW, 1.2, 0 });
		f.push_back(p = { -1 / 2.0 * wallW,0, 0 });

	}

	// floor board
	{
		double
			angle = -60.0 / 360 * 2 * PI,
			R = abs(100 / sin(angle)),
			x0 = -R * cos(angle),
			y0 = 120;
		vector<vector<double>> poly;

		std::function<double(double)> func = [&](double y) {
			y -= y0;
			return sqrt(pow(R + 15, 2) - y * y) + x0;
		};

		{
			poly.push_back({ 0, 0, 0 });
			poly.push_back({ 120, 0, 0 });
			poly.push_back({ 120, 20,0 });
			poly.push_back({ 114, 20,0 });
			poly.push_back({ 114, 16,0 });
			poly.push_back({ 110, 16,0 });
			poly.push_back({ 110, 20,0 });
		}
		{
			int y1 = 20, y2 = 50;

			for (int i = 0; i < 360; i++) {
				double
					y = i / 360.0 * (y2 - y1) + y1,
					x = func(y);

				poly.push_back({ x, y ,0 });
			}
		}
		{
			poly.push_back({ 120, 50, 0 });
			poly.push_back({ 120, 70, 0 });
			poly.push_back({ 114, 70,0 });
			poly.push_back({ 114, 66,0 });
			poly.push_back({ 110, 66,0 });
			poly.push_back({ 110, 70,0 });
		}
		{
			int y1 = 70, y2 = 100;

			for (int i = 0; i < 360; i++) {
				double
					y = i / 360.0 * (y2 - y1) + y1,
					x = func(y);

				poly.push_back({ x, y ,0 });
			}
		}
		{
			poly.push_back({ 120, 100, 0 });
			poly.push_back({ 120, 120, 0 });
			poly.push_back({ 114, 120,0 });
			poly.push_back({ 114, 116,0 });
			poly.push_back({ 110, 116,0 });
			poly.push_back({ 110, 120,0 });
		}

		{
			for (int i = 0; i < 720; i++) {
				double a = i / 720.0 * angle;

				poly.push_back({
					R * cos(a) + x0,
					R * sin(a) + y0, 0 });
			}
		}
		poly.push_back({ 0, 0, 0 });

		for (int i = 0; i <= 8; i++) {
			md.Translator(p = { 0,0,-0.25 / 2 + i * 4.0 }, p1 = { 0,0,0.25 / 2 + i * 4.0 }, poly);
		}

		//md.Translator(poly, f);
	}

	// Corridor
	{
		double w = 6, l = 30, h = 0.25;

		for (int i = 1; i <= 8 - 1; i+= 2) {
			for (int k = 0; k < 2; k++) {
				md.Cuboid(p = { 120 - w / 2, (k == 0 ? 20 : 70) + l / 2, i * 4.0 }, w, l, h);
				// Handrail
				md.Cuboid(p = { 120, (k == 0 ? 20 : 70) + l / 2, i * 4.0 + 1.2 }, 0.08, l, 0.08);
				md.Cuboid(p = { 120 - w, (k == 0 ? 20 : 70) + l / 2, i * 4.0 + 1.2 }, 0.08, l, 0.08);
				md.Cuboid(p = { 120, (k == 0 ? 20 : 70) + l / 2, i * 4.0 + 0.8 }, 0.04, l, 0.04);
				md.Cuboid(p = { 120 - w, (k == 0 ? 20 : 70) + l / 2, i * 4.0 + 0.8 }, 0.04, l, 0.04);
				md.Cuboid(p = { 120, (k == 0 ? 20 : 70) + l / 2, i * 4.0 + 0.4 }, 0.04, l, 0.04);
				md.Cuboid(p = { 120 - w, (k == 0 ? 20 : 70) + l / 2, i * 4.0 + 0.4 }, 0.04, l, 0.04);
			}
		}
	}
	
	// Stair
	{
		double L = 4, W = 4, H = 4, PlatformL = 1;

		Stair_3(md, p = { 120 - 6 - W / 2, 20 - L, 0.25 / 2 }, p1 = {0, 1, 0}, L, W, H, 8, 30, PlatformL);
		Stair_3(md, p = { 120 - 6 - W / 2, 70 - L, 0.25 / 2 }, p1 = {0, 1, 0}, L, W, H, 8, 30, PlatformL);
		Stair_3(md, p = { 120 - 6 - W / 2,120 - L, 0.25 / 2 }, p1 = {0, 1, 0}, L, W, H, 8, 30, PlatformL);
	}
	
}


void Library_Building_Beiyou_1(Modeling& md) {

	// floor board
	{
		vector<vector<double>> poly;
		Modeling::Point p(3), p1(3);
		double boardH = 0.25, H = 5;

		vector<Modeling::Point> f;
		{
			double wallW = 0.2;

			f.push_back(p = { -1 / 2.0 * wallW,0, 0 });
			f.push_back(p = { 1 / 2.0 * wallW,0, 0 });
			f.push_back(p = { 1 / 2.0 * wallW, 1.2, 0 });
			f.push_back(p = { -1 / 2.0 * wallW, 1.2, 0 });
			f.push_back(p = { -1 / 2.0 * wallW,0, 0 });

		}

		{
			// Layer 1st
			poly.push_back({ -90, 0, 0 });
			poly.push_back({ -90, -160, 0 });
			poly.push_back({ -50, -160, 0 });
			poly.push_back({ -50, -80, 0 });
			poly.push_back({ 50, -80, 0 });
			poly.push_back({ 50, -160, 0 });
			poly.push_back({ 90, -160, 0 });
			poly.push_back({ 90, 0, 0 });
			poly.push_back({ -90, 0, 0 });

			md.Translator(p = { 0,0, -boardH / 2 }, p1 = { 0,0, boardH / 2 }, poly);
			md.Translator(p = { 0,0, -boardH / 2 + H }, p1 = { 0,0, boardH / 2 + H }, poly);
			Stair(md, p = { 0, -80 - 20, 0 }, p1 = { 0, 1, 0 }, 20, 50, H, 40);
		}

		{// Layer 2+st
			for (int i = 2; i <= 10; i++) {
				poly.clear();

				double t = std::max(-160.0, -160 + 10.0 * (i - 2));

				poly.push_back({ -90, 0, 0 });
				poly.push_back({ -50, -20, 0 });
				poly.push_back({ -20, -20, 0 });
				poly.push_back({ -20, -60, 0 });
				poly.push_back({ -50, -60, 0 });
				poly.push_back({ -50, -20, 0 });
				poly.push_back({ -90,  0, 0 });

				poly.push_back({ -90, t, 0 });
				poly.push_back({ -50, t, 0 });
				poly.push_back({ -50, -80, 0 });

				if (i == 2) {
					double R = 25;

					for (int i = 0; i < 360; i++) {
						double a = PI - i / 360.0 * PI;
						poly.push_back({ R * cos(a), R * sin(a) - 80, 0 });
					}
				}


				poly.push_back({ 50, -80, 0 });
				poly.push_back({ 50, t, 0 });
				poly.push_back({ 90, t, 0 });
				poly.push_back({ 90, 0, 0 });

				poly.push_back({ 50, -20, 0 });
				poly.push_back({ 50, -60, 0 });
				poly.push_back({ 20, -60, 0 });
				poly.push_back({ 20, -20, 0 });
				poly.push_back({ 50, -20, 0 });
				poly.push_back({ 90,   0, 0 });

				md.Translator(p = { 0,0, -boardH / 2 + i * H }, p1 = { 0,0, boardH / 2 + i * H }, poly);

			}
		}
	}
	// Stair
	{
		double L = 4, W = 4, H = 5, PlatformL = 1;
		Modeling::Point p(3), p1(3);

		Stair_2(md,
			p = { 20 + W / 2, -40 + 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, 1, 0 },
			L, W, H, 9, 40, PlatformL);

		Stair_2(md,
			p = { 20 + W / 2, -40 - 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, -1, 0 },
			L, W, H, 9, 40, PlatformL, -1);

		Stair_2(md,
			p = { -20 - W / 2, -40 + 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, 1, 0 },
			L, W, H, 9, 40, PlatformL, -1);

		Stair_2(md,
			p = { -20 - W / 2, -40 - 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, -1, 0 },
			L, W, H, 9, 40, PlatformL);


		for (int i = 2; i <= 10; i++) {
			md.Cuboid(p = { 20 + W / 2, -40, i * H }, 4, 2, 0.25);
			md.Cuboid(p = { -20 - W / 2, -40, i * H }, 4, 2, 0.25);
		}

		for (int i = 2; i <= 10; i++) {
			double W = 2;

			Stair(md, p = { -50 + W / 2, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);
			Stair(md, p = { -50 + W / 2 * 3, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);
			Stair(md, p = { 50 - W / 2, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);
			Stair(md, p = { 50 - W / 2 * 3, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);

			if (i != 10) md.Cuboid(p = { 50 - W, -40 - 5 - 1, i * H }, 2 * W, 2, 0.25);
			md.Cuboid(p = { 50 - W, -40 + 5 + 1, i * H }, 2 * W, 2, 0.25);
			if (i != 10) md.Cuboid(p = { -50 + W, -40 - 5 - 1, i * H }, 2 * W, 2, 0.25);
			md.Cuboid(p = { -50 + W, -40 + 5 + 1, i * H }, 2 * W, 2, 0.25);
		}

		for (int i = 3; i <= 10; i++) {
			double W = 10, L = 6;
			double t = std::max(-160.0, -160 + 10.0 * (i - 2));

			Stair(md, p = { -90 + W / 2, t - L, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, L, W, H, 40);
			Stair(md, p = { 90 - W / 2, t - L, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, L, W, H, 40);
		}
	}

	// bookshelf
	{
	}
}

void Library_Building_Beiyou_1(Modeling& md, Modeling& md_glass) {

	// floor board
	{
		vector<vector<double>> poly;
		Modeling::Point p(3), p1(3);
		double boardH = 0.25, H = 5;

		vector<Modeling::Point> f;

		{
			// Layer 1st
			poly.push_back({ -90, 0, 0 });
			poly.push_back({ -90, -160, 0 });
			poly.push_back({ -50, -160, 0 });
			poly.push_back({ -50, -80, 0 });
			poly.push_back({  50, -80, 0 });
			poly.push_back({  50, -160, 0 });
			poly.push_back({  90, -160, 0 });
			poly.push_back({  90, 0, 0 });
			poly.push_back({ -90, 0, 0 });

			md.Translator(p = { 0,0, -boardH / 2 }, p1 = { 0,0, boardH / 2 }, poly);
			md.Translator(p = { 0,0, -boardH / 2 + H }, p1 = { 0,0, boardH / 2 + H }, poly);
			Stair(md, p = { 0, -80 - 20, 0 }, p1 = {0, 1, 0}, 20, 50, H, 40);
		}

		{// Layer 2+st
			for (int i = 2; i <= 10; i++) {
				poly.clear();

				double t = std::max(-160.0, -160 + 10.0 * (i - 2));
				
				poly.push_back({ -90, 0, 0 });
				poly.push_back({ -50, -20, 0 });
				poly.push_back({ -20, -20, 0 });
				poly.push_back({ -20, -60, 0 });
				poly.push_back({ -50, -60, 0 });
				poly.push_back({ -50, -20, 0 });
				poly.push_back({ -90,  0, 0 });

				poly.push_back({ -90, t, 0});
				poly.push_back({ -50, t, 0 });
				poly.push_back({ -50, -80, 0 });

				if (i == 2) {
					double R = 25;

					for (int i = 0; i < 360; i++) {
						double a = PI - i / 360.0 * PI;
						poly.push_back({ R * cos(a), R * sin(a) - 80, 0 });
					}
				}


				poly.push_back({ 50, -80, 0 });
				poly.push_back({ 50, t, 0 });
				poly.push_back({ 90, t, 0 });
				poly.push_back({ 90, 0, 0 });

				poly.push_back({ 50, -20, 0 });
				poly.push_back({ 50, -60, 0 });
				poly.push_back({ 20, -60, 0 });
				poly.push_back({ 20, -20, 0 });
				poly.push_back({ 50, -20, 0 });
				poly.push_back({ 90,   0, 0 });
			
				md.Translator(p = { 0,0, -boardH / 2 + i * H }, p1 = { 0,0, boardH / 2 + i * H }, poly);

			}
		}
	}
	// Stair
	{
		double L = 4, W = 4, H = 5, PlatformL = 1;
		Modeling::Point p(3), p1(3);

		Stair_2(md, 
			p = { 20 + W / 2, -40 + 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, 1, 0 }, 
			L, W, H, 9, 40, PlatformL);

		Stair_2(md,
			p = { 20 + W / 2, -40 - 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, -1, 0 },
			L, W, H, 9, 40, PlatformL, -1);

		Stair_2(md,
			p = { -20 - W / 2, -40 + 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, 1, 0 },
			L, W, H, 9, 40, PlatformL, -1);

		Stair_2(md,
			p = { -20 - W / 2, -40 - 2 / 2.0, H + 0.25 / 2 },
			p1 = { 0, -1, 0 },
			L, W, H, 9, 40, PlatformL);


		for (int i = 2; i <= 10; i++) {
			md.Cuboid(p = { 20 + W / 2, -40, i * H }, 4, 2, 0.25);
			md.Cuboid(p = {-20 - W / 2, -40, i * H }, 4, 2, 0.25);
		}

		for (int i = 2; i <= 10; i++) {
			double W = 2;

			Stair(md, p = { -50 + W / 2, -40 - 5, (i -1)  * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);
			Stair(md, p = { -50 + W / 2 * 3, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);
			Stair(md, p = { 50 - W / 2, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);
			Stair(md, p = { 50 - W / 2 * 3, -40 - 5, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, 10, 2, H, 40);

			if (i != 10) md.Cuboid(p = { 50 - W, -40 - 5 - 1, i * H }, 2 * W, 2, 0.25);
			md.Cuboid(p = { 50 - W, -40 + 5 + 1, i * H }, 2 * W, 2, 0.25);
			if (i != 10) md.Cuboid(p = { -50 + W, -40 - 5 - 1, i * H }, 2 * W, 2, 0.25);
			md.Cuboid(p = { -50 + W, -40 + 5 + 1, i * H }, 2 * W, 2, 0.25);
		}

		for (int i = 3; i <= 10; i++) {
			double W = 10, L = 6;
			double t = std::max(-160.0, -160 + 10.0 * (i - 2));

			Stair(md, p = { -90 + W / 2, t - L, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, L, W, H, 40);
			Stair(md, p = {  90 - W / 2, t - L, (i - 1) * H + 0.25 / 2 },
				p1 = { 0, 1, 0 }, L, W, H, 40);
		}
	}

	// bookshelf
	{
	}
}

void Bookshelf(Modeling& md, Modeling::Point certer, double L, double W, double H, double later) {

}

void Library_Building_1(Modeling& md) {
	double R = 90, H = 5, h = 0.25;
	int N = 1800, n = 12;

	vector<Modeling::Point> f;
	Modeling::Point p(3), p1(3), p2(3);

	f.push_back({ R / 2, 0, 0 });
	f.push_back({ R, 0, 0 });
	f.push_back({ R, h, 0 });
	f.push_back({ R / 2, h, 0 });
	f.push_back({ R / 2, 0, 0 });

	md.Rotator_Translator(p = { 0, 0, 0 }, p1 = { 0, 0, 1 }, f, 
		p2 = {0, 0, 1}, n * 5, 720 * n, 0, n * 2 * PI
	);

	{
		f.clear();
		f.push_back({ R / 2 - 6, 0, 0 });
		f.push_back({ R / 2 - 6, 0, 0 });
		f.push_back({ R / 2 - 6, h, 0 });
		f.push_back({ R / 2, h, 0 });
		f.push_back({ R / 2, 0, 0 });

		vector<Modeling::Point> r, r1, r2;
		r.push_back({ R / 2 - 6 - 0.08, 0, 0 });
		r.push_back({ R / 2 - 6, 0, 0 });
		r.push_back({ R / 2 - 6, 0.08, 0 });
		r.push_back({ R / 2 - 6 - 0.08, 0.08, 0 });
		r.push_back({ R / 2 - 6 - 0.08, 0, 0 });

		r1.push_back({ R / 2 - 6 - 0.04, 0, 0 });
		r1.push_back({ R / 2 - 6, 0, 0 });
		r1.push_back({ R / 2 - 6, 0.04, 0 });
		r1.push_back({ R / 2 - 6 - 0.04, 0.04, 0 });
		r1.push_back({ R / 2 - 6 - 0.04, 0, 0 });

		r2.push_back({ R / 2 - 6 - 0.04, 0, 0 });
		r2.push_back({ R / 2 - 6, 0, 0 });
		r2.push_back({ R / 2 - 6, 0.04, 0 });
		r2.push_back({ R / 2 - 6 - 0.04, 0.04, 0 });
		r2.push_back({ R / 2 - 6 - 0.04, 0, 0 });

		for (int i = 0; i < 8; i++) {
			md.Rotator_Translator(p = { 0, 0, 0 }, p1 = { 0, 0, 1 }, f,
				p2 = { 0, 0, 1 }, n * 5, 720 * n, PI / 4 * i, PI / 4 * i + PI
			);
			md.Rotator_Translator(p = { 0, 0, 1.2 }, p1 = { 0, 0, 1 }, r,
				p2 = { 0, 0, 1 }, n * 5, 720 * n, PI / 4 * i, PI / 4 * i + PI
			);
			md.Rotator_Translator(p = { 0, 0, 0.8 }, p1 = { 0, 0, 1 }, r1,
				p2 = { 0, 0, 1 }, n * 5, 720 * n, PI / 4 * i, PI / 4 * i + PI
			);
			md.Rotator_Translator(p = { 0, 0, 0.4 }, p1 = { 0, 0, 1 }, r2,
				p2 = { 0, 0, 1 }, n * 5, 720 * n, PI / 4 * i, PI / 4 * i + PI
			);
		}
	}

}