#ifndef GRAPHICS_PLOT_H
#define GRAPHICS_PLOT_H

#include "../../../../Math/src/Matrix/Mat.h"
#include "RGB.h"

namespace Graphics {

	/* 
	 * 画等高线
	 */
	void contour(Mat<ARGB>& image, Mat<>& map, const int N) {
		int x_step[] = { 1, 0, 1 },
			y_step[] = { 0, 1, 1 };
		double 
			max = map.max(),
			min = map.min(),								//get the max & min of the map
			delta = (max - min) / N, 
			layer = min;

		for (int i = 0; i <= N; i++, layer += delta) {		//for N layer between max & min, get the edge of each layer
			for (int y = 0; y < map.rows - 1; y++) {		//for every point(x,y) to compute
				for (int x = 0; x < map.cols - 1; x++) {
					int flag = map.data[y * map.cols + x] >= layer ? 1 : 0;
					for (int k = 0; k < 3; k++) {			//basic unit is 2x2 matrix
						int xt = x + x_step[k],
							yt = y + y_step[k];
						if (
							(map.data[yt * map.cols + xt] >= layer ? 1 : 0) != flag
							) {
							flag = 2; break;
						}
					}
					if (flag == 2) {
						for (int k = 0; k < 3; k++) {
							int xt = x + x_step[k],
								yt = y + y_step[k];
							if (map.data[yt * map.cols + xt] >= layer) 
								drawPoint(image, xt, yt);
						}
					}
				}
			}
		}
	}

	void contour(Mat<ARGB>& image, Mat<>& map) {
		double 
			min   = map.min(),
			delta = map.max() - min;

		for (int i = 0; i < map.size(); i++)
			image(map.i2x(i), map.i2y(i)) = colorlist((map[i] - min) / delta, 1);
	}

	void contour(Mat<ARGB>& image, Mat<>& mapX, Mat<>& mapY, Mat<>& mapZ) {
		double
			minX = mapX.min(), maxX = mapX.max(),
			minY = mapY.min(), maxY = mapY.max(),
			minZ = mapZ.min(), maxZ = mapZ.max();

		for (int i = 0; i < mapX.size(); i++) {
			ARGB color = 0;
			color += (ARGB)((mapX[i] - minX) / (maxX - minX) * 0xFF) * 0x10000;
			color += (ARGB)((mapY[i] - minY) / (maxY - minY) * 0xFF) * 0x100;
			color += (ARGB)((mapZ[i] - minZ) / (maxZ - minZ) * 0xFF) * 0x1;
			image(mapX.i2x(i), mapX.i2y(i)) = color;
		}
	}

	/*---------------- 色谱 ----------------*/
	ARGB colorlist(double index, int model = 1)
	{
		double A = 0, R = 0, G = 0, B = 0, a = index, b = 1 - a;
		switch (model)
		{
		case 1: {
			B = a <= 9.0 / 16 ? (a < 1.0 / 16 ? 0.5 + 8 * a : (a > 6.0 / 16 ? 1 - (16 / 3.0) * (a - 6.0 / 16) : 1)) : 0;
			R = b <= 9.0 / 16 ? (b < 1.0 / 16 ? 0.5 + 8 * b : (b > 6.0 / 16 ? 1 - (16 / 3.0) * (b - 6.0 / 16) : 1)) : 0;
			G = (a >= 3.0 / 16 && b >= 3.0 / 16) ? (a < 6.0 / 16 ? (16 / 3.0) * (a - 3.0 / 16) : (b < 6.0 / 16 ? (16 / 3.0) * (b - 3.0 / 16) : 1)) : 0;
		}break;
		case 2: {
			B = a <= 9.0 / 16 ? (a < 1.0 / 16 ? 0.5 + 8 * a : (a > 6.0 / 16 ? 1 - (16 / 3.0) * (a - 6.0 / 16) : 1)) : 0;
			R = b <= 9.0 / 16 ? (b < 1.0 / 16 ? 0.5 + 8 * b : (b > 6.0 / 16 ? 1 - (16 / 3.0) * (b - 6.0 / 16) : 1)) : 0;
			G = (a >= 3.0 / 16 && b >= 3.0 / 16) ? (a < 6.0 / 16 ? (16 / 3.0) * (a - 3.0 / 16) : (b < 6.0 / 16 ? (16 / 3.0) * (b - 3.0 / 16) : 1)) : 0;
			A = 0.8;
		}break;
		}
		A *= 0xFF, R *= 0xFF, G *= 0xFF, B *= 0xFF;
		return (ARGB)A * 0x1000000 + (ARGB)R * 0x10000 + (ARGB)G * 0x100 + (ARGB)B;
	}
}

#endif