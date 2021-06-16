/*
Copyright 2020,2021 LiGuer. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/
#include "Plot.h"

/*--------------------------------[ plot ]--------------------------------*/
void Plot::plot(Mat<>& x, Mat<>& y) {
	static bool isinit = true;
	if (isinit) {
		isinit = false;
		pmin.alloc(2) = { x.min(),y.min() }; 
		pmax.alloc(2) = { x.max(),y.max() };
		pdiff.sub(pmax, pmin);
		Mat<> tmp; tmp.mul(0.1, pdiff);
		G.setAxisLim(pmin -= tmp, pmax += tmp); pmin += tmp; pmax -= tmp;
	}
	for (int i = 0; i < x.size() - 1; i++) 
		G.drawLine(
			x[i], x[i + 1],
			y[i], y[i + 1]
		);
}
void Plot::plot(Mat<>& x, Mat<>& y, Mat<>& z) {
	for (int i = 0; i < x.size() - 1; i++)
		G.drawLine(
			x[i], x[i + 1],
			y[i], y[i + 1],
			z[i], z[i + 1]
		);
}
/*--------------------------------[ drawAxis ]--------------------------------*/
void Plot::axis() {
	Mat<> p(2);
	G.drawRectangle(pmin, pmax); 
	G.g.FontSize = 10;
	for (int y = pmin[1]; y <= pmax[1]; y++)
		for (int x = pmin[0]; x <= pmax[0]; x++)
			G.drawLine(x, x, pmin[1], pmin[1] + pdiff[1] / 100),
			G.drawLine(x, x, pmax[1] - pdiff[1] / 100, pmax[1]),
			G.drawLine(pmin[0], pmin[0] + pdiff[0] / 100, y, y),
			G.drawLine(pmax[0] - pdiff[0] / 100, pmax[0], y, y),
			G.drawNum(p = { (double)x, pmin[1] - pdiff[1] / 100 }, x),
			G.drawNum(p = { pmin[0] - pdiff[0] / 20, (double)y }, y);
}
void Plot::grid() {
	G.g.PaintColor = 0xFC000000;
	for (int y = pmin[1]; y <= pmax[1]; y++)
		for (int x = pmin[0]; x <= pmax[0]; x++)
			G.drawLine(x, x, pmin[1], pmax[1]),
			G.drawLine(pmin[0], pmax[0], y, y);
	G.g.PaintColor = 0x0;
}
void Plot::title(const char* words, int n) {
	Mat<> p(2);
	G.drawString(p = { (pmin[0] + pmax[0]) / 2 - pdiff[0] / 100 * n, pmax[1] + pdiff[1] / 20 }, words, n);
}