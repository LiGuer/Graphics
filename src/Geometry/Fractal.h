#ifndef FRACTAL_H
#define FRACTAL_H

#include <math.h>
#include <complex>
#include <vector>
#include "../Matrix/Mat.h"

#define PI 3.141692653589

using namespace std;

namespace Fractal {
/*
 *  Mandelbrot Set
 */
// 是否属于MandelbrotSet/JuliaSet
int isMandelbrotSet(std::complex<double> C, std::complex<double> Z0, int iterateTimes) {
	std::complex<double> z = Z0;
	for (int epoch = 0; epoch < iterateTimes; epoch++) {	// 迭代
		if (abs(z) > 2) return epoch;						// |Zn|>2不可能收敛
		z = z * z + C;										// Zn+1 = Zn² + C   
	} return 0;							//属于,输出0; 不属于,输出判断出的当次迭代数
}

void Mandelbrot(std::complex<double> min, std::complex<double> max, int resSize, int iterateTimes, Mat<int> Set) {
	Set.zero(resSize, resSize);
	double deltaReal = (max.real() - min.real()) / resSize,
	       deltaImag = (max.imag() - min.imag()) / resSize;
	for (int i = 0; i < resSize; i++) {
		for (int j = 0; j < resSize; j++) {
			double real = min.real() + deltaReal * i,
				   imag = min.imag() + deltaImag * j;
			Set(i, j) = isMandelbrotSet(
				std::complex<double>(real, imag), 
				std::complex<double>(0, 0), 
				iterateTimes
			);
		}
	}
}

/*
 *  Julia集
 */
void Julia(std::complex<double> C, std::complex<double> min, std::complex<double> max, int resSize, int iterateTimes, Mat<int> Set) {
	Set.zero(resSize, resSize);
	double deltaReal = (max.real() - min.real()) / resSize,
	       deltaImag = (max.imag() - min.imag()) / resSize;
	for (int i = 0; i < resSize; i++) {
		for (int j = 0; j < resSize; j++) {
			double real = min.real() + deltaReal * i, 
				   imag = min.imag() + deltaImag * j;
			Set(i, j) = isMandelbrotSet(C, std::complex<double>(real, imag), iterateTimes);
		}
	}
}


/*
 *  Hilbert 曲线
 */
int HilbertCurve_xy2d(int n, int x, int y) {	//convert (x,y) to d
	int distance = 0;
	for (int s = n / 2; s > 0; s /= 2) {
		int rx = (x & s) > 0, 
			ry = (y & s) > 0;
		distance += s * s * ((0b11 * rx) ^ ry);
		// rotation
		if (ry == 0) {
			if (rx == 1) { 
				x = n - 1 - x; 
				y = n - 1 - y;
			}
			int t = x; x = y; y = t;
		}
	}return distance;
}
void HilbertCurve_d2xy(int n, int distance, int& x, int& y) {//convert d to (x,y)
	x = y = 0;
	for (int s = 1; s < n; s *= 2) {
		int rx = 1 & (distance / 2), 
			ry = 1 & (distance ^ rx);
		// rotation
		if (ry == 0) {
			if (rx == 1) { 
				x = n - 1 - x; 
				y = n - 1 - y;
			}
			int t = x; x = y; y = t;
		}
		x += s * rx; 
		y += s * ry;
		distance /= 4;
	}
}


/*
 *  Perlin Noise
 */
double PerlinNoise(double x, double y, Mat<>& randomGridGradient) {
	// 对四个格点
	int x0[] = { x, x + 1, x, x + 1 }, 
		y0[] = { y, y, y + 1, y + 1 };
	double n[4];
	for (int i = 0; i < 4; i++) {
		//[1][2]
		double random = randomGridGradient(x0[i], y0[i]);
		n[i] = (x - x0[i]) * cos(random) 
			 + (y - y0[i]) * sin(random);
	}
	//[3]
	double 
		sx = x - (int)x,
		sy = y - (int)y,
		ix0 = (n[1] - n[0]) * (3.0 - sx * 2.0) * sx * sx + n[0],
		ix1 = (n[3] - n[2]) * (3.0 - sx * 2.0) * sx * sx + n[2];
	return (ix1 - ix0) * (3.0 - sy * 2.0) * sy * sy + ix0;
}
Mat<>& PerlinNoise(Mat<>& output, int frequency) {
	Mat<> randomGridGradient;
	randomGridGradient.rands(frequency + 1, frequency + 1, 0, 256);
	for (int y = 0; y < output.cols; y++)
		for (int x = 0; x < output.rows; x++)
			output(x, y) = PerlinNoise(
				(double)x / output.rows * frequency,
				(double)y / output.cols * frequency,
				randomGridGradient
			);
	return output;
}

}
#endif