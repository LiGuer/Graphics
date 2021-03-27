/*
Copyright 2020 LiGuer. All Rights Reserved.
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
#include "Graphics3D.h"
Mat<double> Graphics3D::TransformMat;											//变换矩阵
/******************************************************************************

*                    底层函数

******************************************************************************/
/*--------------------------------[ 初始化 ]--------------------------------*/
void Graphics3D::init(int width, int height) {
	g.init(width, height);
	Z_index.zero(height, width);
	for (int i = 0; i < Z_index.rows * Z_index.cols; i++)Z_index.data[i] = -0x7FFFFFFF;
	TransformMat.E(4);
}
/*--------------------------------[ value2pix ]--------------------------------*/
void Graphics3D::value2pix(int x0, int y0, int z0, int& x, int& y, int& z) {
	Mat<double> point(4, 1);
	{double t[] = { x0, y0, z0, 1 }; point.getData(t); }
	point.mult(TransformMat, point);
	x = point[0]; y = point[1]; z = point[2];
}
void Graphics3D::value2pix(Mat<double>& p0, int& x, int& y, int& z) {
	value2pix(p0[0], p0[1], p0.rows < 3 ? 0 : p0[2], x, y, z);
}
/******************************************************************************

*                    Draw

******************************************************************************/
/*--------------------------------[ 画点 ]--------------------------------*/
void Graphics3D::drawPoint(double x0, double y0, double z0) {
	int x, y, z;
	value2pix(x0, y0, z0, x, y, z);
	if (g.judgeOutRange(y, x) || z < Z_index(y, x))return;
	g.drawPoint(y, x); Z_index(y, x) = z;
}
void Graphics3D::drawPoint(Mat<double>& p0) {
	int x, y, z;
	value2pix(p0, x, y, z);
	if (g.judgeOutRange(y, x) || z < Z_index(y, x))return;
	g.drawPoint(y, x); Z_index(y, x) = z;
}
/*--------------------------------[ 画直线 ]--------------------------------*/
void Graphics3D::drawLine(double sx0, double ex0, double sy0, double ey0, double sz0, double ez0) {
	int sx, sy, sz, ex, ey, ez;
	value2pix(sx0, sy0, sz0, sx, sy, sz); value2pix(ex0, ey0, ez0, ex, ey, ez);
	int err[3] = { 0 }, inc[3] = { 0 };
	int delta[3] = { ex - sx, ey - sy, ez - sz };						//计算坐标增量
	int index[3] = { sx, sy, sz };
	//设置xyz单步方向	
	for (int dim = 0; dim < 3; dim++) {
		if (delta[dim] > 0)inc[dim] = 1; 								//向右
		else if (delta[dim] == 0)inc[dim] = 0;							//垂直
		else { inc[dim] = -1; delta[dim] = -delta[dim]; }				//向左
	}
	int distance = delta[0] > delta[1] ? delta[0] : delta[1];//总步数
	distance = delta[2] > distance ? delta[2] : distance;//总步数
	//画线
	for (int i = 0; i <= distance + 1; i++) {
		//唯一输出：画点

		if (!g.judgeOutRange(index[1], index[0]) && index[2] >= Z_index(index[1], index[0])) {
			g.drawPoint(index[1], index[0]); Z_index(index[1], index[0]) = index[2];
		}
		for (int dim = 0; dim < 3; dim++) {						//xyz走一步
			err[dim] += delta[dim];
			if (err[dim] > distance) { err[dim] -= distance; index[dim] += inc[dim]; }
		}
	}
}
void Graphics3D::drawLine(Mat<double>& sp, Mat<double>& ep) {
	int sx, sy, sz, ex, ey, ez;
	value2pix(sp, sx, sy, sz); value2pix(ep, ex, ey, ez);

	int err[3] = { 0 }, inc[3] = { 0 };
	int delta[3] = { ex - sx, ey - sy, ez - sz };						//计算坐标增量
	int index[3] = { sx, sy, sz };
	//设置xyz单步方向	
	for (int dim = 0; dim < 3; dim++) {
		if (delta[dim] > 0)inc[dim] = 1; 								//向右
		else if (delta[dim] == 0)inc[dim] = 0;							//垂直
		else { inc[dim] = -1; delta[dim] = -delta[dim]; }				//向左
	}
	int distance = delta[0] > delta[1] ? delta[0] : delta[1];//总步数
	distance = delta[2] > distance ? delta[2] : distance;//总步数
	//画线
	for (int i = 0; i <= distance + 1; i++) {
		//唯一输出：画点

		if (!g.judgeOutRange(index[1], index[0]) && index[2] >= Z_index(index[1], index[0])){ 
			g.drawPoint(index[1], index[0]); Z_index(index[1], index[0]) = index[2]; 
		}
		for (int dim = 0; dim < 3; dim++) {						//xyz走一步
			err[dim] += delta[dim];
			if (err[dim] > distance) { err[dim] -= distance; index[dim] += inc[dim]; }
		}
	}
}
/*--------------------------------[ 画折线 ]--------------------------------*/
void Graphics3D::drawPolyline(Mat<double> *p, int n) {
	for (int i = 0; i < n - 1; i++) drawLine(p[i], p[i + 1]);
}
/*--------------------------------[ 画等高线 ]--------------------------------*/
void Graphics3D::contour(Mat<double>& map, const int N) {
	int x_step[] = { 1,0,1 }, y_step[] = { 0,1,1 };
	double max = map.max(), min = map.min();			//get the max & min of the map
	double delta = (max - min) / N, layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {		//for N layer between max & min, get the edge of each layer
		for (int y = 0; y < map.rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map.cols - 1; x++) {
				int flag = map.data[y * map.cols + x] >= layer ? 1 : 0;
				for (int k = 0; k < 3; k++) {			//basic unit is 2x2 matrix
					int xt = x + x_step[k], yt = y + y_step[k];
					int flagtemp = map.data[yt * map.cols + xt] >= layer ? 1 : 0;
					if (flagtemp != flag) { flag = 2; break; }
				}
				if (flag == 2) {
					for (int k = 0; k < 3; k++) {
						int xt = x + x_step[k], yt = y + y_step[k];
						if (map.data[yt * map.cols + xt] >= layer) drawPoint(xt, yt);
					}
				}
			}
		}
	}
}
void Graphics3D::contourface(Mat<double>& map, const int N) {
	double max = map.max(), min = map.min();			//get the max & min of the map
	double delta = (max - min) / N, layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {		//for N layer between max & min, get the edge of each layer
		for (int y = 0; y < map.rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map.cols - 1; x++) {
				if (map.data[y * map.cols + x] >= layer)
					g.setPoint(x, y, colorlist((double)i/N, 1));
			}
		}
	}
}
/*--------------------------------[ 画三角形 ]--------------------------------*/
void Graphics3D::drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3) {
	drawLine(p1, p2); drawLine(p2, p3); drawLine(p3, p1);
}
/*--------------------------------[ 画矩形 ]--------------------------------*/
void Graphics3D::drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct) {
	if (direct == NULL) {
		Mat<double> pt(2, 1);
		{ double t[] = { sp[0],ep[1] }; pt.getData(t); }
		drawLine(sp, pt); drawLine(pt, ep);
		{ double t[] = { sp[1],ep[0] }; pt.getData(t); }
		drawLine(sp, pt); drawLine(pt, ep);
	}
}
void Graphics3D::drawRectangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4) {
	drawLine(p1, p2); drawLine(p2, p3); drawLine(p3, p4); drawLine(p4, p1);
}
/*--------------------------------[ 画圆 ]--------------------------------
*	[约束方程]:
		平面点法式: A(x-x0) + B(y-y0) + C(z-z0) = 0    n=(A,B,C)
		球方程: (x-x0)² + (y-y0)² + (z-z0)² = r²
			x = r cosΦ cosθ + x0    Φ∈[-π/2, π/2]
			y = r cosΦ sinθ + y0    θ∈[-π, π]
			z = r sinΦ + z0
	[推导]:
		A cosΦ cosθ + B cosΦ sinθ + C sinΦ = 0
		对于某一θ值:
		(A cosθ + B sinθ)cosΦ  + C sinΦ = 0
		C1 cosΦ + C sinΦ = 0
		sin(Φ + α) = 0    α = arcsin(C1 / sqrt(C1² + C²))
		Φ = - arcsin(C1 / sqrt(C1² + C²))
**-----------------------------------------------------------------------*/
void Graphics3D::drawCircle(Mat<double>& center, double r, Mat<double>* direct) {
}
/*--------------------------------[ 画椭圆 ]--------------------------------
*	[约束方程]:
		平面点法式: A(x-x0) + B(y-y0) + C(z-z0) = 0    n=(A,B,C)
		球方程: (x-x0)² + (y-y0)² + (z-z0)² = r²
			x = r cosΦ cosθ + x0    Φ∈[-π/2, π/2]
			y = r cosΦ sinθ + y0    θ∈[-π, π]
			z = r sinΦ + z0
	[推导]:
		A cosΦ cosθ + B cosΦ sinθ + C sinΦ = 0
		对于某一θ值:
		(A cosθ + B sinθ)cosΦ  + C sinΦ = 0
		C1 cosΦ + C sinΦ = 0
		sin(Φ + α) = 0    α = arcsin(C1 / sqrt(C1² + C²))
		Φ = - arcsin(C1 / sqrt(C1² + C²))
**-----------------------------------------------------------------------*/
void Graphics3D::drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct) {
}
/*--------------------------------[ 画多边形 ]--------------------------------*/
void Graphics3D::drawPolygon(Mat<double> p[],int n) {
	for (int i = 0; i < n; i++) drawLine(p[i], p[(i + 1) % n]);
}
/*--------------------------------[ 画曲面 ]--------------------------------*/
void Graphics3D::drawSurface(Mat<double> z, double xs, double xe, double ys, double ye) {
	Mat<double> p(3, 1), pl(3, 1), pu(3, 1);
	double dx = (xe - xs) / z.rows, dy = (ye - ys) / z.cols;
	for (int y = 0; y < z.cols; y++) {
		for (int x = 0; x < z.rows; x++) {
			{double t[] = { xs + x * dx,ys + y * dy,z(x,y) }; p.getData(t); }
			if (x > 0) { 
				double t[] = { xs + (x - 1) * dx,ys + y * dy,z(x - 1,y) }; 
				pl.getData(t); drawLine(pl, p);
			}
			if (y > 0) { 
				double t[] = { xs + x * dx,ys + (y - 1) * dy,z(x,y - 1) }; 
				pu.getData(t); drawLine(pu, p);
			}
		}
	}
}
/*--------------------------------[ 画四面体 ]--------------------------------*/
void Graphics3D::drawTetrahedron(Mat<double> p[]) {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < i; j++)
			drawLine(p[i], p[j]);
}
/*--------------------------------[ 画矩体 ]--------------------------------
*	矩体共十二条边，从对角点引出6条: 
		(x0,y0,z0)&(x1,y0,z0)  (x0,y0,z0)&(x0,y1,z0)  (x0,y0,z0)&(x0,y0,z1)
		(x1,y1,z1)&(x0,y1,z1)  (x1,y1,z1)&(x1,y0,z1)  (x1,y1,z1)&(x1,y1,z0)
	另外六条:
		(x1,y0,z0)&(x0,y1,z0)  (x1,y0,z0)&(x0,y0,z1)
		(x0,y1,z0)&(x0,y0,z1)  (x0,y1,z0)&(x1,y0,z0)
		(x0,y0,z1)&(x1,y0,z0)  (x0,y0,z1)&(x0,y1,z1)
**------------------------------------------------------------------------*/
void Graphics3D::drawCuboid(Mat<double> pMin, Mat<double> pMax) {
	Mat<double> pMinTemp, pMaxTemp;
	for (int i = 0; i < 3; i++) {
		pMinTemp = pMin; pMaxTemp = pMax;
		pMinTemp[i] = pMax[i]; pMaxTemp[(i + 1) % 3] = pMin[(i + 1) % 3];
		drawLine(pMin, pMinTemp);
		drawLine(pMax, pMaxTemp);
		drawLine(pMinTemp, pMaxTemp);
		pMaxTemp[(i + 1) % 3] = pMax[(i + 1) % 3]; pMaxTemp[(i + 2) % 3] = pMin[(i + 2) % 3];
		drawLine(pMinTemp, pMaxTemp);
	}
}
/*--------------------------------[ 画圆台 ]--------------------------------
**------------------------------------------------------------------------*/
void Graphics3D::drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta) {
	// 计算 Rotate Matrix
	Mat<double> direction, rotateAxis, rotateMat(4), tmp(3, 1);
	direction.add(st, ed.negative(direction));
	{double t[] = { 0, 0, 1 }; tmp.getData(t); }
	rotateAxis.crossProduct(direction, tmp);
	double rotateAngle = -acos(tmp.dot(direction, tmp) / direction.norm());
	{
		rotate(rotateAxis, rotateAngle, tmp.zero(3, 1), rotateMat);
		tmp = rotateMat; rotateMat.E(3);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				rotateMat(i, j) = tmp(i, j);
	}
	// 画圆台
	Mat<double> stPoint, edPoint, preStPoint, preEdPoint, deltaVector(3, 1);
	for (int i = 0; i <= 360 / delta; i++) {
		{double t[] = { cos(i * delta * 2.0 * PI / 360), sin(i * delta * 2.0 * PI / 360),0 }; deltaVector.getData(t); }
		deltaVector.mult(rotateMat, deltaVector);
		stPoint.add(st, stPoint.mult(Rst, deltaVector));
		edPoint.add(ed, edPoint.mult(Red, deltaVector));
		if (i != 0) {
			drawLine(stPoint, preStPoint); drawLine(st, stPoint);
			drawLine(edPoint, preEdPoint); drawLine(ed, edPoint);
			drawLine(stPoint, edPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}
/*--------------------------------[ 画圆柱 ]--------------------------------*/
void Graphics3D::drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta) {
	drawFrustum(st, ed, r, r, delta);
}
/*--------------------------------[ 画球 ]--------------------------------
*	[公式]: x² + y² + z² = R²
		参数方程,点集:
			x = r cosΦ cosθ    Φ∈[-π/2, π/2]
			y = r cosΦ sinθ    θ∈[0, 2π]
			z = r sinΦ
*	[流程]:
		[1] 画纬度线
		[2] 画经度线
**-----------------------------------------------------------------------*/
void Graphics3D::drawSphere(Mat<double>& center, double r, int delta) {
	// 经纬度法
	Mat<double> point(3, 1);
	for (int i = 0; i < 360 / delta; i++) {
		double theta = (i * delta) * 2.0 * PI / 360;
		for (int j = -90 / delta; j <= 90 / delta; j++) {
			double phi = (j * delta) * 2.0 * PI / 360;
			point[0] = r * cos(phi) * cos(theta) + center[0];
			point[1] = r * cos(phi) * sin(theta) + center[1];
			point[2] = r * sin(phi) + center[2];
			drawPoint(point);
		}
	}
}
/*--------------------------------[ getSphereFibonacciPoint 球面均匀点分布 ]--------------------------------
*	[Referance]:
		[1] Thanks and copyright for https://github.com/SebLague/Boids
**---------------------------------------------------------------------------------------------------------*/
void Graphics3D::drawSphere2(Mat<double>& center, double r, int n) {
	// 均匀球面点
	Mat<double> point(3, 1);
	double goldenRatio = (1 + sqrt(5)) / 2;				// 黄金分割点
	double angleIncrement = PI * 2 * goldenRatio;
	for (int i = 0; i < 300; i++) {
		double t = (double)i / n, inclination = acos(1 - 2 * t), azimuth = angleIncrement * i;
		point[0] = center[0] + r * sin(inclination) * cos(azimuth);
		point[1] = center[1] + r * sin(inclination) * sin(azimuth);
		point[2] = center[2] + r * cos(inclination);
		drawPoint(point);
	}
}
/*--------------------------------[ 画椭球 ]--------------------------------
*	[公式]: (x/rx)² + (y/ry)² + (z/rz)² = 1
		参数方程,点集:
			x = rx cosΦ cosθ    Φ∈[-π/2, π/2]
			y = ry cosΦ sinθ    θ∈[-π, π]
			z = rz sinΦ
*	[流程]:
		[1] 画纬度线
		[2] 画经度线
**-----------------------------------------------------------------------*/
void Graphics3D::drawEllipsoid(Mat<double>& center, Mat<double>& r) {
	const int delta = 5;
	Mat<double> point(3, 1);
	for (int i = 0; i < 360 / delta; i++) {
		double theta = (i * delta) * 2.0 * PI / 360;
		for (int j = -90 / delta; j <= 90 / delta; j++) {
			double phi = (j * delta) * 2.0 * PI / 360;
			point[0] = r[0] * cos(phi) * cos(theta) + center[0];
			point[1] = r[1] * cos(phi) * sin(theta) + center[1];
			point[2] = r[2] * sin(phi) + center[2];
			drawPoint(point);
		}
	}
}
/*--------------------------------[ 网格 ]--------------------------------*/
void Graphics3D::grid() {
	//int x0 = value2pix(0, 0), y0 = value2pix(0, 1);		//原点像素坐标
	//double size[2] = { pSizeMax[0] - pSizeMin[0],pSizeMax[1] - pSizeMin[1] };
	int x0, y0;
	double size[2];
	/*------ 网格 ------*/
	g.PaintColor = 0x00ccff;
	g.PaintSize = 1;
	/*------ 计算间隔值 ------*/
	size[0] /= 10; size[1] /= 10;
	double delta[2] = { 1,1 };
	for (int dim = 0; dim < 2; dim++) {
		while ((int)size[dim] == 0) {
			size[dim] *= 10;
			delta[dim] /= 10;
		}
		while ((int)size[dim] >= 10) {
			size[dim] /= 10;
			delta[dim] *= 10;
		}
		delta[dim] *= (int)size[dim];
	}
	//g.drawGrid(x0, y0, 0, 0, -value2pix(delta[0], 0), -value2pix(delta[1], 1));
	//g.drawGrid(x0, y0, g.gWidth, 0, value2pix(delta[0], 0), -value2pix(delta[1], 1));
	//g.drawGrid(x0, y0, 0, g.gHeight, -value2pix(delta[0], 0), value2pix(delta[1], 1));
	//g.drawGrid(x0, y0, g.gWidth, g.gHeight, value2pix(delta[0], 0), value2pix(delta[1], 1));
	/*------ 坐标轴 ------*/
	g.PaintColor = 0xffffff;
	g.PaintSize = 3;
	//g.drawLine(x0, coor2pix(pSizeMin[1], 1), x0, coor2pix(pSizeMax[1], 1));
	//g.drawLine(coor2pix(pSizeMin[0], 0), y0, coor2pix(pSizeMax[0], 0), y0);
	/*------ 轴标号 ------*/
	g.PaintSize = 1;
	g.FontSize = 50;
	g.PaintColor = 0xffffff;
	g.drawChar(x0 - 40, y0 + 15, '0');
}
/*---------------- 色谱 ----------------*/
Graphics::ARGB Graphics3D::colorlist(double index, int model)
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
	return (Graphics::ARGB)A * 0x1000000 + (Graphics::ARGB)R * 0x10000 + (Graphics::ARGB)G * 0x100 + (Graphics::ARGB)B;
}

/******************************************************************************

*                    Transformation - 3D

******************************************************************************/
/*--------------------------------[ 平移 ]--------------------------------
[x'] = [ 1  0  0  dx] [x]
|y'|   | 0  1  0  dy| |y|
|z'|   | 0  0  1  dz| |z|
[1 ]   [ 0  0  0  1 ] [1]
**-----------------------------------------------------------------------*/
void Graphics3D::translation(Mat<double>& delta, Mat<double>& transMat) {
	int n = delta.rows;
	Mat<double> translationMat(n + 1);
	for (int i = 0; i < n; i++)translationMat(i, n) = delta[i];
	transMat.mult(translationMat, transMat);
}
/*--------------------------------[ 三维旋转·四元数 ]--------------------------------
*	[公式]: v' = q v q`¹
		q = [cos(θ/2), u s in(θ/2)]
		v=>[0,v]经旋转轴u旋转Ѳ角后得到v'
	多次旋转:
		v' = q1q2 v q2`¹q1`¹ = (q1q2) v (q1q2)`¹
	四元数化旋转矩阵:
		四元数左乘:
		q v =	[a -b -c -d] v
				|b  a -d  c|
				|c  d  a -b|
				[d -c  b  a]
		四元数右乘:
		v q =	[a -b -c -d] v
				|b  a  d -c|
				|c -d  a  b|
				[d  c -b  a]
**--------------------------------------------------------------------------*/
void Graphics3D::rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& transMat) {
	Mat<double> tmp;
	translation(center.negative(tmp), transMat);
	// rotate
	rotateAxis.normalization();
	Mat<double> q(4, 1);
	q[0] = cos(theta / 2);
	for (int i = 1; i < 4; i++) q[i] = rotateAxis[i - 1] * sin(theta / 2);
	Mat<double> rotateMat(4, 4);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			rotateMat(i, j) = q[((j % 2 == 0 ? 1 : -1) * i + j + 4) % 4];
	for (int i = 1; i < 4; i++)rotateMat(i, i % 3 + 1) = -rotateMat(i, i % 3 + 1);
	Mat<double> rotateMatR(rotateMat);
	for (int i = 1; i < 4; i++) {
		rotateMat(0, i) = -rotateMat(0, i); rotateMatR(i, 0) = -rotateMatR(i, 0);
	}
	rotateMat.mult(rotateMat, rotateMatR);
	rotateMatR = rotateMat;
	rotateMat.E(4);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			rotateMat(i, j) = rotateMatR(i + 1, j + 1);
	transMat.mult(rotateMat, transMat);
	translation(center, transMat);
}
/*--------------------------------[ 缩放 ]--------------------------------
[x'] = [ sx 0  0  0 ] [x]
|y'|   | 0  sy 0  0 | |y|
|z'|   | 0  0  sz 0 | |z|
[1 ]   [ 0  0  0  1 ] [1]
**-----------------------------------------------------------------------*/
void Graphics3D::scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat) {
	Mat<double> tmp;
	translation(center.negative(tmp), transMat);
	// scaling
	int n = scale.rows;
	Mat<double> scaleMat(n + 1);
	for (int i = 0; i < n; i++)scaleMat(i, i) = scale[i];
	transMat.mult(scaleMat, transMat);
	translation(center, transMat);
}
/******************************************************************************

*                    Set

******************************************************************************/
void Graphics3D::setAxisLim(Mat<double> pMin, Mat<double> pMax) {
	Mat<double> center, scale, zero(pMin.rows, 1);
	center.add(pMax, pMin.negative(center));
	scale = center;
	center.mult(1.0 / 2, center);
	translation(center);
	for (int i = 0; i < scale.rows; i++) {
		if (i == 1)scale[i] = g.Canvas.rows / scale[i];
		else scale[i] = g.Canvas.cols / scale[i];
	}
	scaling(scale, center);
}