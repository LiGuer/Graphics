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
/******************************************************************************
*                    Transformation-3D
******************************************************************************/
/*--------------------------------[ 平移 ]--------------------------------
[x'] = [ 1  0  0  dx] [x]
|y'|   | 0  1  0  dy| |y|
 z'|   | 0  0  1  dz| |z|
[1 ]   [ 0  0  0  1 ] [1]
------------------------------------------------------------------------*/
void Graphics3D::translation(Mat<double>& delta) {
	int n = delta.rows;
	Mat<double> translationMat(n + 1);
	for (int i = 0; i < n; i++)translationMat(i, n) = delta[i];
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
void Graphics3D::rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center) {
	Mat<double> temp;
	translation(center.negative(temp));
	// rotate
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
	TransformMat.mult(rotateMat, TransformMat);

	translation(center);
}
/*--------------------------------[ 缩放 ]--------------------------------
[x'] = [ sx 0  0  0 ] [x]
|y'|   | 0  sy 0  0 | |y|
 z'|   | 0  0  sz 0 | |z|
[1 ]   [ 0  0  0  1 ] [1]
**-----------------------------------------------------------------------*/
void Graphics3D::scaling(Mat<double>& scale, Mat<double>& center) {
	Mat<double> temp;
	translation(center.negative(temp));
	// scaling
	int n = scale.rows;
	Mat<double> scalingMat(n + 1);
	for (int i = 0; i < n; i++)scalingMat(i, i) = scale[i];
	TransformMat.mult(scalingMat, TransformMat);

	translation(center);
}