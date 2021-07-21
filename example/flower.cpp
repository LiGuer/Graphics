#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"

int main() {
	GraphicsND G(1000, 1000);

	double r, o = 0, R = 250;
	std::vector<double> os, as,as2; Mat<> t(3), zero(3); zero[2] = R;
	for (int i = 0; i < 90; i++) {
		o += (PI / 2) / 90;
		r = 3 * 100 * sin(o) * cos(o) / (pow(sin(o), 3) + pow(cos(o), 3));
		os.push_back(o);
		as.push_back(asin(r / R) - PI / 2);
		as2.push_back(asin(r / (2*R)) - PI / 2);
	}
	
	{
		G.FaceColor = 0xFFAEC8; G.isLineTriangleSet = 1;
		zero[2] = R;
		for (int i = 0; i < os.size(); i++)
			for (int j = 0; j < 6; j++)
				G.drawSphere(zero, R, os[i] + PI * 2 / 6 * j, os[i] + PI * 2 / 6 * j + PI / 2 / 90 * 2, -PI / 2, as[i], PI / 2 / 90);
		zero[2] = R/3;
		for (int i = 0; i < os.size(); i++)
			for (int j = 0; j < 6; j++)
				G.drawSphere(zero, R/3, os[i] + PI * 2 / 12 + PI * 2 / 6 * j, os[i] + PI * 2 / 12 + PI * 2 / 6 * j + PI / 2 / 90 * 2, -PI / 2, 2*(as[i] + PI/2) - PI/2, PI / 2 / 90);
		zero[2] = 2*R;
		for (int i = 0; i < os.size(); i++)
			for (int j = 0; j < 6; j++)
				G.drawSphere(zero, 2*R, os[i] + PI * 2 / 12 + PI * 2 / 6 * j, os[i] + PI * 2 / 12 + PI * 2 / 6 * j + PI / 2 / 90 * 2, -PI / 2, as2[i], PI / 2 / 90);
		G.writeModel("D:/LiGu.stl");
	}
}