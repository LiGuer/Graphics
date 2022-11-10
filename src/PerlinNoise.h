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