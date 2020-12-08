﻿/*
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
/******************************************************************************
*                    底层函数
******************************************************************************/
/*--------------------------------[ 构造函数 ]--------------------------------*/
Graphics3D::Graphics3D(int WindowSize_Width, int WindowSize_height) {
	init(WindowSize_Width, WindowSize_height);
}
/*--------------------------------[ 析构函数 ]--------------------------------*/
Graphics3D::~Graphics3D() {
	if (g != NULL)free(g);
}
/*--------------------------------[ 初始化 ]--------------------------------*/
void Graphics3D::init(int WindowSize_Width, int WindowSize_height) {
	if (g != NULL)free(g);
	g = new Graphics;
	g->init(WindowSize_Width, WindowSize_height);
	TransformMat.E(4);
}
/******************************************************************************
*                    Draw
******************************************************************************/
/*--------------------------------[ 画点 ]--------------------------------*/
void Graphics3D::drawPoint(Mat<double>& p0) {
	Mat<double> point(4, 1);
	for (int i = 0; i < 3; i++)point[i] = p0[i]; point[3] = 1;
	point.mult(TransformMat, point);
	g->drawPoint((int)point[0], (int)point[1]);
}
/*--------------------------------[ 画直线 ]--------------------------------*/
void Graphics3D::drawLine(Mat<double>& sp, Mat<double>& ep) {
	Mat<double> startpoint(4, 1), endpoint(4, 1);
	for (int i = 0; i < 3; i++) {
		startpoint[i] = sp[i];
		endpoint[i] = ep[i];
	}startpoint[3] = 1; endpoint[3] = 1;
	startpoint.mult(TransformMat, startpoint);
	endpoint.mult(TransformMat, endpoint);
	g->drawLine(startpoint[0], startpoint[1], endpoint[0], endpoint[1]);
}
/*--------------------------------[ 画多边形 ]--------------------------------*/
void Graphics3D::drawPolygon(Mat<double> p[],int n) {
	for (int i = 0; i < n; i++) drawLine(p[i], p[(i + 1) % n]);
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
void Graphics3D::drawCircle(Mat<double>& center, double r, Mat<double>& direction) {
	Mat<double> point[72];
	for (int i = 0; i < 72; i++) {
		point[i].zero(3, 1);
		double theta = (i * 5.0) * 2.0 * PI / 360;
		double C1 = direction[0] * cos(theta) + direction[1] * sin(theta);
		double Ct = sqrt(C1 * C1 + direction[2] * direction[2]);
		double phi = -asin(C1 / Ct);
		point[i][0] = r * cos(phi) * cos(theta) + center[0];
		point[i][1] = r * cos(phi) * sin(theta) + center[1];
		point[i][2] = r * sin(phi) + center[2];
	}
	drawPolygon(point, 72);
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
void Graphics3D::drawSphere(Mat<double>& center, double r) {
	// 经纬度法
	const int delta = 5;
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
/*--------------------------------[ 三角填充 ]--------------------------------
**-----------------------------------------------------------------------*/
void Graphics3D::fillTriangle(Mat<double> p0[]) {
	int x[3], y[3];
	Mat<double> point(4, 1);
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < 3; i++)point[i] = p0[k][i]; point[3] = 1;
		point.mult(TransformMat, point);
		x[k] = point[0]; y[k] = point[1];
	}
	g->fillPolygon(x, y, 3);
}
/******************************************************************************
*                    Transformation-3D
******************************************************************************/
/*--------------------------------[ 平移 ]--------------------------------
[x'] = [ 1  0  0  dx] [x]
|y'|   | 0  1  0  dy| |y|
|z'|   | 0  0  1  dz| |z|
[1 ]   [ 0  0  0  1 ] [1]
**-----------------------------------------------------------------------*/
void Graphics3D::translation(Mat<double>& delta, Mat<double>& translationMat) {
	int n = delta.rows;
	translationMat.E(n + 1);
	for (int i = 0; i < n; i++)translationMat(i, n) = delta[i];
}
void Graphics3D::translation(Mat<double>& delta) {
	Mat<double> translationMat;
	translation(delta, translationMat);
	TransformMat.mult(translationMat, TransformMat);
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
void Graphics3D::rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& rotateMat) {
	Mat<double> temp, translationMat;
	translation(center.negative(temp), translationMat);
	// rotate
	Mat<double> q(4, 1);
	q[0] = cos(theta / 2);
	for (int i = 1; i < 4; i++) q[i] = rotateAxis[i - 1] * sin(theta / 2);
	rotateMat.zero(4, 4);
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

	rotateMat.mult(rotateMat, translationMat);
	translation(center, translationMat);
	rotateMat.mult(rotateMat, translationMat);
}
void Graphics3D::rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center) {
	Mat<double> rotateMat;
	rotate(rotateAxis, theta, center, rotateMat);
	TransformMat.mult(rotateMat, TransformMat);
}
/*--------------------------------[ 缩放 ]--------------------------------
[x'] = [ sx 0  0  0 ] [x]
|y'|   | 0  sy 0  0 | |y|
|z'|   | 0  0  sz 0 | |z|
[1 ]   [ 0  0  0  1 ] [1]
**-----------------------------------------------------------------------*/
void Graphics3D::scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& scaleMat) {
	Mat<double> temp, translationMat;
	translation(center.negative(temp), translationMat);
	// scaling
	int n = scale.rows;
	scaleMat.E(n + 1);
	for (int i = 0; i < n; i++)scaleMat(i, i) = scale[i];

	scaleMat.mult(scaleMat, translationMat);
	translation(center, translationMat);
	scaleMat.mult(scaleMat, translationMat);
}
void Graphics3D::scaling(Mat<double>& scale, Mat<double>& center) {
	Mat<double> scaleMat;
	scaling(scale, center, scaleMat);
	TransformMat.mult(scaleMat, TransformMat);
}