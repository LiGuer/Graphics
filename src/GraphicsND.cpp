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
#include "GraphicsND.h"
Mat<double> GraphicsND::TransformMat;											//变换矩阵
/******************************************************************************

*                    底层函数

******************************************************************************/
/*--------------------------------[ 初始化 ]--------------------------------*/
void GraphicsND::init(int width, int height) {
	g.init(width, height);
	Z_Buffer.zero(height, width);
	for (int i = 0; i < Z_Buffer.rows * Z_Buffer.cols; i++)Z_Buffer.data[i] = -0x7FFFFFFF;
	TransformMat.E(4);
}
/*--------------------------------[ value2pix ]--------------------------------*/
void GraphicsND::value2pix(int x0, int y0, int z0, int& x, int& y, int& z) {
	Mat<double> point(4, 1);
	{double t[] = { x0, y0, z0, 1 }; point.getData(t); }
	point.mult(TransformMat, point);
	x = point[0]; y = point[1]; z = point[2];
}
void GraphicsND::value2pix(Mat<double>& p0, Mat<int>& pAns) {
	pAns.zero(p0.rows, 1);
	Mat<double> point(4, 1);
	{double t[] = { p0[0], p0[1], p0.rows < 3 ? 0 : p0[2], 1 }; point.getData(t); }
	point.mult(TransformMat, point);
	for (int i = 0; i < pAns.rows; i++)pAns[i] = point[i];
}
/*--------------------------------[ setPix ]--------------------------------*/
bool GraphicsND::setPix(int x,int y, int z) {
		if (g.judgeOutRange(y, x) || z < Z_Buffer(y, x))return false;
		g.drawPoint(y, x); Z_Buffer(y, x) = z; return true;
}
/******************************************************************************

*                    Draw

******************************************************************************/
/*--------------------------------[ 画点 ]--------------------------------*/
void GraphicsND::drawPoint(double x0, double y0, double z0) {
	int x, y, z;
	value2pix(x0, y0, z0, x, y, z);
	setPix(x, y, z);
}
void GraphicsND::drawPoint(Mat<double>& p0) {
	Mat<int> p;
	value2pix(p0, p);
	setPix(p[0], p[1], p[2]);
}
/*--------------------------------[ 画直线 ]--------------------------------*/
void GraphicsND::drawLine(double sx0, double ex0, double sy0, double ey0, double sz0, double ez0) {
	int sx, sy, sz, ex, ey, ez;
	value2pix(sx0, sy0, sz0, sx, sy, sz); value2pix(ex0, ey0, ez0, ex, ey, ez);
	int err[3] = { 0 }, inc[3] = { 0 }, delta[3] = { ex - sx, ey - sy, ez - sz }, index[3] = { sx, sy, sz };
	//设置xyz单步方向	
	for (int dim = 0; dim < 3; dim++) {
		inc[dim] = delta[dim] == 0 ? 0 : (delta[dim] > 0 ? 1 : -1);		//符号函数(向右,垂直,向左)
		delta[dim] *= delta[dim] < 0 ? -1 : 1;							//向左
	}
	int distance = delta[0] > delta[1] ? delta[0] : delta[1];//总步数
	distance = delta[2] > distance ? delta[2] : distance;//总步数
	//画线
	for (int i = 0; i <= distance + 1; i++) {
		setPix(index[0], index[1], index[2]);					//唯一输出：画点
		for (int dim = 0; dim < 3; dim++) {						//xyz走一步
			err[dim] += delta[dim];
			if (err[dim] > distance) { err[dim] -= distance; index[dim] += inc[dim]; }
		}
	}
}
void GraphicsND::drawLine(Mat<double>& sp0, Mat<double>& ep0) {
	Mat<int> sp, ep;
	value2pix(sp0, sp); value2pix(ep0, ep);
	Mat<int> err(sp.rows, 1), inc(sp.rows, 1), delta, index(sp), tmp;
	delta.add(ep, sp.negative(tmp));
	//设置xyz单步方向	
	for (int dim = 0; dim < sp.rows; dim++) {
		inc[dim] = delta[dim] == 0 ? 0 : (delta[dim] > 0 ? 1 : -1);		//符号函数(向右,垂直,向左)
		delta[dim] *= delta[dim] < 0 ? -1 : 1;							//向左
	} int distance = delta.max();										//总步数
	//画线
	for (int i = 0; i <= distance + 1; i++) {
		setPix(index[0], index[1], index[2]);							//唯一输出：画点
		for (int dim = 0; dim < sp.rows; dim++) {						//xyz走一步
			err[dim] += delta[dim];
			if (err[dim] > distance) { err[dim] -= distance; index[dim] += inc[dim]; }
		}
	}
}
/*--------------------------------[ 画折线 ]--------------------------------*/
void GraphicsND::drawPolyline(Mat<double> *p, int n, bool close) {
	for (int i = 0; i < n - 1; i++) drawLine(p[i], p[i + 1]);
	if (close) drawLine(p[0], p[n - 1]);
}
/*--------------------------------[ 画三角形 ]--------------------------------*/
void GraphicsND::drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, bool FACE, bool LINE) {
	Mat<int> pt[3], tmp;
	value2pix(p1, pt[0]); value2pix(p2, pt[1]); value2pix(p3, pt[2]);
	int Dim = p1.rows;
	//调整顺序
	int pUindex = 1;
	for (int i = 1; i < 3; i++) pUindex = pt[i][0] < pt[pUindex][0] ? i : pUindex;
	tmp = pt[0]; pt[0] = pt[pUindex]; pt[pUindex] = tmp;
	//设置xyz单步方向	
	Mat<int> err[2], inc[2], delta[2], index[2];
	int distance[2];
	for (int k = 0; k < 2; k++) {
		err[k].zero(Dim, 1), inc[k].zero(Dim, 1); delta[k].add(pt[k + 1], pt[0].negative(tmp)); index[k] = pt[0];
		for (int dim = 0; dim < Dim; dim++) {
			inc[k][dim] = delta[k][dim] == 0 ? 0 : (delta[k][dim] > 0 ? 1 : -1);//符号函数(向右,垂直,向左)
			delta[k][dim] *= delta[k][dim] < 0 ? -1 : 1;						//向左
		} distance[k] = delta[k].max();											//总步数
	}
	//画线
	int distanceMaxCur = distance[0] > distance[1] ? 0 : 1;
	for (int i = 0; i <= distance[distanceMaxCur] + 1; i++) {
		if (FACE) {
			g.PaintColor = FaceColor;
			Mat<int> errTmp(Dim, 1), incTmp(Dim, 1), deltaTmp, indexTmp(index[0]);
			deltaTmp.add(index[1], index[0].negative(tmp));
			for (int dim = 0; dim < Dim; dim++) {								//设置xyz单步方向	
				incTmp[dim] = deltaTmp[dim] == 0 ? 0 : (deltaTmp[dim] > 0 ? 1 : -1);//符号函数(向右,垂直,向左)
				deltaTmp[dim] *= deltaTmp[dim] < 0 ? -1 : 1;					//向左
			}
			int distanceTmp = deltaTmp.max();									//总步数
			for (int i = 0; i <= distanceTmp + 1; i++) {						//画线
				setPix(indexTmp[0], indexTmp[1], indexTmp[2]);					//唯一输出：画点
				//setPix(indexTmp[0] + 1, indexTmp[1], indexTmp[2]);				//唯一输出：画点
				for (int dim = 0; dim < Dim; dim++) {							//xyz走一步
					errTmp[dim] += deltaTmp[dim];
					if (errTmp[dim] > distanceTmp) { errTmp[dim] -= distanceTmp; indexTmp[dim] += incTmp[dim]; }
				}
			}
		}
		if (LINE) {
			g.PaintColor = LineColor;
			setPix(index[0][0], index[0][1], index[0][2]);
			setPix(index[1][0], index[1][1], index[1][2]);
		}
		for (int dim = 0; dim < Dim; dim++) {						//xyz走一步
			err[0][dim] += delta[0][dim];
			err[1][dim] += delta[1][dim];
			if (err[0][dim] > distance[0]) { err[0][dim] -= distance[0]; index[0][dim] += inc[0][dim]; }
			if (err[1][dim] > distance[1]) { err[1][dim] -= distance[1]; index[1][dim] += inc[1][dim]; }
		}
		// 到达节点
		if (i == distance[1 - distanceMaxCur]) {
			int kt = 1 - distanceMaxCur;
			delta[kt].add(pt[distanceMaxCur + 1], pt[kt + 1].negative(tmp));
			index[kt] = pt[kt + 1];
			for (int dim = 0; dim < Dim; dim++) {
				inc[kt][dim] = delta[kt][dim] == 0 ? 0 : (delta[kt][dim] > 0 ? 1 : -1);	//符号函数(向右,垂直,向左)
				delta[kt][dim] *= delta[kt][dim] < 0 ? -1 : 1;					//向左
			}distance[kt] = delta[kt].max();											//总步数
		}
	}
}
/*--------------------------------[ 画矩形 ]--------------------------------*/
void GraphicsND::drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct) {
	if (direct == NULL) {
		Mat<double> pt(2, 1);
		{ double t[] = { sp[0],ep[1] }; pt.getData(t); }
		drawLine(sp, pt); drawLine(pt, ep);
		{ double t[] = { sp[1],ep[0] }; pt.getData(t); }
		drawLine(sp, pt); drawLine(pt, ep);
	}
	else {
		// 计算 Rotate Matrix
	}
}
void GraphicsND::drawRectangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4) {
	drawLine(p1, p2); drawLine(p2, p3); drawLine(p3, p4); drawLine(p4, p1);
}
/*--------------------------------[ 画多边形 ]--------------------------------*/
void GraphicsND::drawPolygon(Mat<double> p[], int n) { drawPolyline(p, n, true); }
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
void GraphicsND::drawCircle(Mat<double>& center, double r, Mat<double>* direct) {
	// 计算 Rotate Matrix
	//画圆
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
void GraphicsND::drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct) {
}
/*--------------------------------[ 画曲面 ]--------------------------------*/
void GraphicsND::drawSurface(Mat<double> z, double xs, double xe, double ys, double ye) {
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
/*--------------------------------[ 画等高线 ]--------------------------------*/
void GraphicsND::contour(Mat<double>& map, const int N) {
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
void GraphicsND::contourface(Mat<double>& map, const int N) {
	double max = map.max(), min = map.min();			//get the max & min of the map
	double delta = (max - min) / N, layer = min;
	for (int i = 0; i <= N; i++, layer += delta) {		//for N layer between max & min, get the edge of each layer
		for (int y = 0; y < map.rows - 1; y++) {		//for every point(x,y) to compute
			for (int x = 0; x < map.cols - 1; x++) {
				if (map.data[y * map.cols + x] >= layer)
					g.setPoint(x, y, colorlist((double)i / N, 1));
			}
		}
	}
}
/*--------------------------------[ 画四面体 ]--------------------------------*/
void GraphicsND::drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4) {
	drawLine(p1, p2); drawLine(p1, p3); drawLine(p1, p4);
	drawLine(p2, p3); drawLine(p2, p4);
	drawLine(p3, p4);
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
void GraphicsND::drawCuboid(Mat<double>& pMin, Mat<double>& pMax) {
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
void GraphicsND::drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta, bool FACE, bool LINE) {
	// 计算 Rotate Matrix
	Mat<double> direction, rotateAxis, rotateMat(4), zAxis(3, 1), tmp; {double t[] = { 0, 0, 1 }; zAxis.getData(t); }
	direction.add(ed, st.negative(direction));
	if (direction[0] != 0 || direction[1] != 0) {
		rotate(
			rotateAxis.crossProduct(direction, zAxis),
			-acos(tmp.dot(direction, zAxis) / direction.norm()),
			tmp.zero(3, 1), rotateMat
		);
		rotateMat.cut(0, 2, 0, 2, rotateMat);
	} else rotateMat.E(3);
	// 画圆台
	Mat<double> stPoint, edPoint, preStPoint, preEdPoint, deltaVector(3, 1);
	for (int i = 0; i <= 360 / delta; i++) {
		{double t[] = { cos(i * delta * 2.0 * PI / 360), sin(i * delta * 2.0 * PI / 360),0 }; deltaVector.getData(t); }
		deltaVector.mult(rotateMat, deltaVector);
		stPoint.add(st, stPoint.mult(Rst, deltaVector));
		edPoint.add(ed, edPoint.mult(Red, deltaVector));
		if (i != 0) {
			if (FACE) {
				g.PaintColor = FaceColor;
				drawTriangle(stPoint, preStPoint, edPoint, true, false);
				drawTriangle(preStPoint, preEdPoint, edPoint, true, false);
			}
			if (LINE) {
				g.PaintColor = LineColor;
				drawLine(stPoint, preStPoint); drawLine(st, stPoint);
				drawLine(edPoint, preEdPoint); drawLine(ed, edPoint);
				drawLine(stPoint, edPoint);
			}
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}
/*--------------------------------[ 画圆柱 ]--------------------------------*/
void GraphicsND::drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta) {
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
void GraphicsND::drawSphere(Mat<double>& center, double r, int delta) {
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
void GraphicsND::drawSphere2(Mat<double>& center, double r, int n) {
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
void GraphicsND::drawEllipsoid(Mat<double>& center, Mat<double>& r) {
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
/*--------------------------------[ 画坐标轴 ]--------------------------------*/
void GraphicsND::drawAxis(double Xmax,double Ymax,double Zmax, bool negative) {
	// 轴
	drawLine(negative ? -Xmax : 0, Xmax);//x
	drawLine(0, 0, negative ? -Ymax : 0, Ymax);//y
	drawLine(0, 0, 0, 0, negative ? -Zmax : 0, Zmax);//z
	// 箭头
	Mat<double> st(3, 1), ed(3, 1);
	int vectorLength = 10, vectorWidth = vectorLength / 2.718281828456;
	{double t[] = { Xmax,0, 0 }; st.getData(t); }
	{double t[] = { Xmax + vectorLength,0, 0 }; ed.getData(t); }
	if (Xmax != 0) drawFrustum(st, ed, vectorWidth, 0, 45);
	{double t[] = { 0,Ymax, 0 }; st.getData(t); }
	{double t[] = { 0,Ymax + vectorLength, 0 }; ed.getData(t); }
	if (Ymax != 0) drawFrustum(st, ed, vectorWidth, 0, 45);
	{double t[] = { 0,0,Zmax }; st.getData(t); }
	{double t[] = { 0,0,Zmax + vectorLength }; ed.getData(t); }
	if (Zmax != 0) drawFrustum(st, ed, vectorWidth, 0, 45);
}
/*--------------------------------[ 网格 ]--------------------------------*/
void GraphicsND::drawGrid() {
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
Graphics::ARGB GraphicsND::colorlist(double index, int model)
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
void GraphicsND::translation(Mat<double>& delta, Mat<double>& transMat) {
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
void GraphicsND::rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& transMat) {
	Mat<double> tmp;
	translation(center.negative(tmp), transMat);
	// rotate
	tmp.mult(sin(theta / 2), rotateAxis.normalization());
	Mat<double> q(4, 1), rotateMat(4); 
	{ double t[] = { cos(theta / 2), tmp[0], tmp[1], tmp[2] }; q.getData(t); }
	for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) rotateMat(i, j) = q[((j % 2 == 0 ? 1 : -1) * i + j + 4) % 4];
	for (int i = 1; i < 4; i++) rotateMat(i, i % 3 + 1) = -rotateMat(i, i % 3 + 1);
	Mat<double> rotateMatR(rotateMat);
	for (int i = 1; i < 4; i++) { rotateMat(0, i) *= -1; rotateMatR(i, 0) *= -1; }
	rotateMat.mult(rotateMat, rotateMatR);
	tmp = rotateMat; rotateMat.E(4);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) rotateMat(i, j) = tmp(i + 1, j + 1);
	transMat.mult(rotateMat, transMat);
	translation(center, transMat);
}
/*--------------------------------[ 缩放 ]--------------------------------
[x'] = [ sx 0  0  0 ] [x]
|y'|   | 0  sy 0  0 | |y|
|z'|   | 0  0  sz 0 | |z|
[1 ]   [ 0  0  0  1 ] [1]
**-----------------------------------------------------------------------*/
void GraphicsND::scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat) {
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
void GraphicsND::setAxisLim(Mat<double> pMin, Mat<double> pMax) {
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