#include "../Matrix/Matrix.h"

void Tree3D(vector<double>& center, 
	        vector<double>& direction, 
	        double length, double angle, double scale,
			int forks, int level,
	        vector<vector<double>>& stSet,
	        vector<vector<double>>& edSet
) {
	if (level < 0)
		return;
	
	Mat<double> rotateMat(3, 3);
	vector<double> new_x(3), new_y(3), p(3), new_direction(3), new_center(3);

	if (1 - abs(dot(direction, p = { 0, 0, 1 })) < 10e-4) 
		new_y = { 0, 1, 0 };
	else 
		normalize(cross(new_y, direction, p = { 0, 0, 1 }));

	normalize(cross(new_x, direction, new_y));
	assign(rotateMat, { new_x, new_y , direction });


	for (int i = 0; i < forks; i++) {
		double theta = 2 * PI * (double)i / forks;
		p = { 
			sin(angle) * cos(theta), 
			sin(angle) * sin(theta),
			cos(angle)
		};

		mul(new_direction, rotateMat, p);
		mul(new_center, length, new_direction);
		add(new_center, center, new_center);

		edSet.push_back(new_center);
		stSet.push_back(center);

		Tree3D(
			new_center, 
			new_direction, 
			length * scale, angle, scale, 
			forks, level - 1, 
			stSet, edSet
		);
	}
}
