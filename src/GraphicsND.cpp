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
void GraphicsND::init(int width, int height, int Dim) { 
	g.init(width, height);  
	Z_Buffer.zero(Dim - 2, 1);
	for (int i = 0; i < Z_Buffer.rows; i++) 
		Z_Buffer[i].zero(g.Canvas.rows, g.Canvas.cols);
	clear(0);
}
void GraphicsND::clear(ARGB color) {
	g.clear(color);
	for (int i = 0; i < Z_Buffer.rows; i++)
		for (int j = 0; j < Z_Buffer[i].rows * Z_Buffer[i].cols; j++)
			Z_Buffer[i].data[j] = -0x7FFFFFFF;
	TransformMat.E(Z_Buffer.rows + 2 + 1);
}
/*--------------------------------[ 点 To 像素 ]--------------------------------*/
void GraphicsND::value2pix(int x0, int y0, int z0, int& x, int& y, int& z) {
	Mat<double> point(TransformMat.rows, 1);
	point[0] = 1; point[1] = x0; point[2] = y0; point[3] = z0;
	point.mult(TransformMat, point);
	x = point[1]; y = point[2]; z = point[3];
}
void GraphicsND::value2pix(Mat<double>& p0, Mat<int>& pAns) {
	pAns.zero(p0.rows, 1);
	Mat<double> point(TransformMat.rows, 1);
	point[0] = 1; for (int i = 0; i < p0.rows; i++) point[i + 1] = p0[i];
	point.mult(TransformMat, point);
	for (int i = 0; i < pAns.rows; i++) pAns[i] = point[i + 1];
}
/*--------------------------------[ 写像素 (正投影) ]--------------------------------*/
bool GraphicsND::setPix(int x,int y, int z, int size) {
	if (g.judgeOutRange(y, x) || z < Z_Buffer[0](y, x))return false;
	if (size == -1) g.drawPoint(y, x);
	else if (size == 0) g.setPoint(y, x, FaceColor); 
	Z_Buffer[0](y, x) = z; return true;
}
bool GraphicsND::setPix(Mat<int>& p0, int size) {
	if (g.judgeOutRange(p0[1], p0[0])) return false;
	for (int i = 2; i < p0.rows; i++)
		if (p0[i] < Z_Buffer[i - 2](p0[1], p0[0])) 
			return false;
	if (size == -1)  g.drawPoint(p0[1], p0[0]); 
	else if (size == 0)  g.setPoint(p0[1], p0[0], FaceColor);
	for (int i = 2; i < p0.rows; i++) Z_Buffer[i - 2](p0[1], p0[0]) = p0[i];
	return true;
}
/*--------------------------------[ 设置坐标范围 ]--------------------------------*/
void GraphicsND::setAxisLim(Mat<double>& pMin, Mat<double>& pMax) {
	Mat<double> scale, tmp;
	scale.add(pMax, pMin.negative(scale));
	for (int i = 0; i < scale.rows; i++) {
		if (i == 1)scale[i] = g.Canvas.rows / scale[i];
		else scale[i] = g.Canvas.cols / scale[i];
	}
	translation(pMin.negative(tmp));
	scaling(scale, tmp.zero(pMin.rows, 1));
}
/******************************************************************************

*                    Draw

******************************************************************************/
/*--------------------------------[ 画点 ]--------------------------------
*	[过程]:
		[1] 原坐标转换至世界参考系的像素坐标
		[2] 绘制像素
------------------------------------------------------------------------*/
void GraphicsND::drawPoint(double x0, double y0, double z0) {
	int x, y, z;
	value2pix(x0, y0, z0, x, y, z);
	setPix(x, y, z);
}
void GraphicsND::drawPoint(Mat<double>& p0) {
	Mat<int> p;
	value2pix(p0, p);
	setPix(p);
}
/*--------------------------------[ 画直线 ]--------------------------------
*	[算法]: Bresenham
*	[特点]:
		1. 化[浮点运算]为[整数运算]：err[]
		2. 各方向均可绘制
*	[原理]:
		设维度边长最大方向 M  (因为设维度长最大方向为扫描方向，才能保证线像素全部填充)
		Δx = (x_max - x_min) / (M_max - M_min) = x_delta / M_delta
		x = x + Δx 
		∵ 离散化只有整数
		若 x = x + 1, 则 1 = n·Δx = n·x_delta / M_delta => n·x_delta = M_delta
		即 x_err 加 x_delta， 直至大于等于 M_delta, 此时 x = x + 1
*	[过程]:
		[1] 原坐标转换至世界参考系的像素坐标
		[2] 计算每个维度的 x_delta, 找到维度长最大方向 M_delta
			通过x_inc 解决(向右,垂直,向左)的方向性问题
		[3] 以 M 维度, 从M_min 至 M_max, 逐像素扫描 (M_i++)
			[4] 画点
			[5] 累计 x_err, 进位 x
---------------------------------------------------------------------------*/
void GraphicsND::drawLine(double sx0, double ex0, double sy0, double ey0, double sz0, double ez0) {
	//[1]
	int sx, sy, sz, ex, ey, ez;
	value2pix(sx0, sy0, sz0, sx, sy, sz); value2pix(ex0, ey0, ez0, ex, ey, ez);
	int err[3] = { 0 }, inc[3] = { 0 }, delta[3] = { ex - sx, ey - sy, ez - sz }, point[3] = { sx, sy, sz };
	//[2]
	for (int dim = 0; dim < 3; dim++) {
		inc[dim] = delta[dim] == 0 ? 0 : (delta[dim] > 0 ? 1 : -1);		//符号函数(向右,垂直,向左)
		delta[dim] *= delta[dim] < 0 ? -1 : 1;							//向左
	}
	int distance = delta[0] > delta[1] ? delta[0] : delta[1];			//总步数
	distance = delta[2] > distance ? delta[2] : distance;				//总步数
	//[3]
	for (int i = 0; i <= distance; i++) {
		setPix(point[0], point[1], point[2]);							//唯一输出：画点
		for (int dim = 0; dim < 3; dim++) {								//xyz走一步
			err[dim] += delta[dim];
			if (err[dim] >= distance) { err[dim] -= distance; point[dim] += inc[dim]; }
		}
	}
}
void GraphicsND::drawLine(Mat<double>& sp0, Mat<double>& ep0) {
	Mat<int> sp, ep;
	value2pix(sp0, sp); value2pix(ep0, ep);
	Mat<int> err(sp.rows, 1), inc(sp.rows, 1), delta, point(sp), tmp;
	delta.add(ep, sp.negative(tmp));
	//设置xyz单步方向	
	for (int dim = 0; dim < sp.rows; dim++) {
		inc[dim] = delta[dim] == 0 ? 0 : (delta[dim] > 0 ? 1 : -1);		//符号函数(向右1,垂直0,向左-1)
		delta[dim] *= delta[dim] < 0 ? -1 : 1;							//绝对值
	} int distance = delta.max();										//总步数
	//画线
	for (int i = 0; i <= distance; i++) {
		setPix(point);													//唯一输出：画点
		for (int dim = 0; dim < sp.rows; dim++) {						//xyz走一步
			err[dim] += delta[dim];
			if (err[dim] >= distance) { err[dim] -= distance; point[dim] += inc[dim]; }
		}
	}
}
/*--------------------------------[ 画折线 ]--------------------------------*/
void GraphicsND::drawPolyline(Mat<double> *p, int n, bool close) {
	for (int i = 0; i < n - 1; i++) drawLine(p[i], p[i + 1]);
	if (close) drawLine(p[0], p[n - 1]);
}
void GraphicsND::drawPolyline(Mat<double>& y, double xmin, double xmax) {
	Mat<double> point(2, 1), prePoint(2, 1);
	{double t[] = { xmin, y[0] }; point.getData(t); prePoint.getData(t); }
	double delta = (xmax - xmin) / (y.cols - 1);
	for (int i = 1; i < y.cols; i++) {
		point[0] += delta; point[1] = y[i];
		drawLine(prePoint, point);
		prePoint = point;
	}
}
/*--------------------------------[ 画三角形 ]--------------------------------
*	[算法]: Bresenham
*	[原理]:
		从 x 最小的顶点 p_x_min 开始,
		以 x 维度为基准, 逐渐跟踪以 p_x_min 为顶点两条边，绘制该两点连接的线
		(绘制线 x 坐标实在相同, 从而避免走样)
		若到达三角形第三顶点，则短边转换至该点方程
*	[过程]:
		[1] 原坐标转换至世界参考系的像素坐标
		[2] 查找 x 最小的顶点 p_x_min
		[3] 类似直线算法, 以 x 方向为基准, 跟踪以 p_x_min 为顶点两条边
			[4] 若到达三角形第三顶点，则短边转换至该点方程
			[5] 画线
-----------------------------------------------------------------------------*/
void GraphicsND::drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, bool FACE, bool LINE) {
	if (FACE) {
		//[1]
		Mat<int> pt[3], tmp;
		value2pix(p1, pt[0]); value2pix(p2, pt[1]); value2pix(p3, pt[2]);
		int Dim = p1.rows;
		//[2]
		int pXminIndex = 0;
		for (int i = 1; i < 3; i++) pXminIndex = pt[i][0] < pt[pXminIndex][0] ? i : pXminIndex;
		if (pXminIndex != 0) { tmp = pt[0]; pt[0] = pt[pXminIndex]; pt[pXminIndex] = tmp; }
		//[3]
		Mat<int> err[2], inc[2], delta[2], point[2];
		for (int k = 0; k < 2; k++) {
			err[k].zero(Dim, 1), inc[k].zero(Dim, 1); delta[k].add(pt[k + 1], pt[0].negative(tmp)); point[k] = pt[0];
			for (int dim = 0; dim < Dim; dim++) {
				inc[k][dim] = delta[k][dim] == 0 ? 0 : (delta[k][dim] > 0 ? 1 : -1);//符号函数(向右,垂直,向左)
				delta[k][dim] *= delta[k][dim] < 0 ? -1 : 1;						//向左
			}
		}
		int dXMaxCur = delta[0][0] > delta[1][0] ? 0 : 1;
		bool flag = true;															//三角形转折点检测开关(不然会二次检测)
		for (int i = 0; i <= delta[dXMaxCur][0]; i++) {
			//[4]三角形转折点检测
			if (i == delta[1 - dXMaxCur][0] && flag) {
				flag = false;														//关闭检测
				int kt = 1 - dXMaxCur;
				delta[kt].add(pt[dXMaxCur + 1], pt[kt + 1].negative(tmp));
				if (delta[kt][0] == 0) break;
				point[kt] = pt[kt + 1];
				for (int dim = 1; dim < Dim; dim++) {
					inc[kt][dim] = delta[kt][dim] == 0 ? 0 : (delta[kt][dim] > 0 ? 1 : -1);	//符号函数(向右,垂直,向左)
					delta[kt][dim] *= delta[kt][dim] < 0 ? -1 : 1;					//向左
				}
			}
			//[5]画线
			Mat<int> errTmp(Dim, 1), incTmp(Dim, 1), deltaTmp, pointTmp(point[0]);
			deltaTmp.add(point[1], point[0].negative(tmp));
			for (int dim = 0; dim < Dim; dim++) {									//设置xyz单步方向	
				incTmp[dim] = deltaTmp[dim] == 0 ? 0 : (deltaTmp[dim] > 0 ? 1 : -1);//符号函数(向右,垂直,向左)
				deltaTmp[dim] *= deltaTmp[dim] < 0 ? -1 : 1;						//向左
			}
			int distanceTmp = deltaTmp.max();										//总步数
			for (int i = 0; i <= distanceTmp; i++) {								//画线
				pointTmp[2]--; setPix(pointTmp, 0);	pointTmp[2]++;					//唯一输出：画点 (Z-1:反走样)
				for (int dim = 0; dim < Dim; dim++) {								//xyz走一步
					errTmp[dim] += deltaTmp[dim];
					if (errTmp[dim] >= distanceTmp) { errTmp[dim] -= distanceTmp; pointTmp[dim] += incTmp[dim]; }
				}
			}
			for (int k = 0; k < 2; k++) {
				point[k][0]++;
				for (int dim = 1; dim < Dim; dim++) {								//xyz走一步
					err[k][dim] += delta[k][dim];
					if (err[k][dim] >= delta[k][0]) { point[k][dim] += err[k][dim] / delta[k][0] * inc[k][dim]; err[k][dim] = err[k][dim] % delta[k][0]; }
				}
			}
		}
	}
	if (LINE) { drawLine(p1, p2); drawLine(p2, p3); drawLine(p3, p1); }
}
/*--------------------------------[ 画矩形 ]--------------------------------*/
void GraphicsND::drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct, bool FACE, bool LINE) {
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
void GraphicsND::drawQuadrilateral(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, bool FACE, bool LINE) {
	if (FACE) {
		drawTriangle(p1, p2, p3, true, false); drawTriangle(p1, p3, p4, FACE, false);
	}
	if (LINE) {
		drawLine(p1, p2); drawLine(p2, p3); drawLine(p3, p4); drawLine(p4, p1);
	}
}
/*--------------------------------[ 画多边形 ]--------------------------------
* [算法]:
		绘制轮数k = 上取整(边数 / 3)
		三角形绘制顶点: i, i + k, i + 2k
* [例子]:
	偶数边:
	四边形: [1] 1 2 3, 3 4 1
	六边形: [1] 1 2 3, 3 4 5, 5 6 1 [2] 1 3 5
	奇数边:
	五边形:	[1] 1 2 3, 3 4 5 [2] 1 3 5
		  
-----------------------------------------------------------------------------*/
void GraphicsND::drawPolygon(Mat<double> p[], int n, bool FACE, bool LINE) { 
	if (FACE) {
		for (int k = 1; k <= (n + 2) / 3; k++) {
			for (int i = 0; i <= n - 2 * k; i += 2 * k) {
				drawTriangle(p[i], p[i + k], p[(i + 2 * k) % n], FACE, false);
			}
		}
	}
	if (LINE) drawPolyline(p, n, true);
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
void GraphicsND::drawCircle(Mat<double>& center, double r, Mat<double>* direct, bool FACE, bool LINE) {
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
/*--------------------------------[ 画四面体 ]--------------------------------*/
void GraphicsND::drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, bool FACE, bool LINE) {
	if (FACE) {
		drawTriangle(p1, p2, p3); drawTriangle(p2, p3, p4); drawTriangle(p3, p4, p1); drawTriangle(p4, p1, p2);
	}
	if (LINE) {
		drawLine(p1, p2); drawLine(p1, p3); drawLine(p1, p4);
		drawLine(p2, p3); drawLine(p2, p4);
		drawLine(p3, p4);
	}
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
void GraphicsND::drawCuboid(Mat<double>& pMin, Mat<double>& pMax, bool FACE, bool LINE) {
	Mat<double> pMinTmp[3], pMaxTmp[3];
	for (int i = 0; i < 3; i++) {
		pMinTmp[i] = pMin; pMinTmp[i][i] = pMax[i];
		pMaxTmp[i] = pMax; pMaxTmp[i][i] = pMin[i];
	}
	if (FACE) {
		for (int i = 0; i < 3; i++) {
			drawQuadrilateral(pMin, pMinTmp[i], pMaxTmp[(i + 2) % 3], pMinTmp[(i + 1) % 3], true, false);
			drawQuadrilateral(pMax, pMaxTmp[i], pMinTmp[(i + 2) % 3], pMaxTmp[(i + 1) % 3], true, false);
		}
	}
	if (LINE) {
		for (int i = 0; i < 3; i++) {
			drawLine(pMin, pMinTmp[i]); drawLine(pMax, pMaxTmp[i]);
			drawLine(pMinTmp[i], pMaxTmp[(i + 1) % 3]); drawLine(pMinTmp[i], pMaxTmp[(i + 2) % 3]);
		}
	}
}
/*--------------------------------[ 画圆台 ]--------------------------------
* [过程]:
		[1] 计算旋转矩阵
		[2] 根据旋转矩阵, 计算绘制点坐标, 完成绘制
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
		rotateMat.cut(1, 3, 1, 3, rotateMat);
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
				drawTriangle(stPoint, preStPoint, edPoint, true, false);
				drawTriangle(preStPoint, preEdPoint, edPoint, true, false);
			}
			if (LINE) {
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
void GraphicsND::drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta, bool FACE, bool LINE) {
	drawFrustum(st, ed, r, r, delta, FACE, LINE);
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
void GraphicsND::drawSphere(Mat<double>& center, double r, int delta, bool FACE, bool LINE) {
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
			if (FACE) {

			}
			if (LINE) {

			}
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
/*--------------------------------[ 画线 any-D ]--------------------------------*/
void GraphicsND::drawSuperLine(Mat<double>* p0, bool FACE, bool LINE) {

}
/*--------------------------------[ 画立方体 any-D ]--------------------------------
*	[算法]: 利用二进制表示各顶点的坐标位置，
			最小顶点为全0，最大顶点为全1
			0: 该维度坐标值 == 最小顶点对应值
			1: 该维度坐标值 == 最大顶点对应值
*	[流程]:
		[1] 以二进制顺序遍历所有顶点
			[2] 连接该点和所有比该点编码多1的点
---------------------------------------------------------------------------*/
void GraphicsND::drawSuperCuboid(Mat<double>& pMin, Mat<double>& pMax, bool FACE, bool LINE) {
	unsigned int Dim = pMin.rows, maxCode = 0;
	Mat<double> st, ed;
	for (int i = 0; i < Dim; i++) maxCode += 1 << i;
	for (unsigned int code = 0; code < maxCode; code++) {
		st = pMin;
		for (int i = 0; i < Dim; i++)
			if (code & (1 << i))
				st[i] = pMax[i];
		ed = st;
		for (int i = 0; i < Dim; i++) {
			if (ed[i] == pMin[i]) {
				ed[i] = pMax[i];
				drawLine(st, ed);
				ed[i] = pMin[i];
			}
		}
	}
}
/*--------------------------------[ 画球 any-D ]--------------------------------
*	[定义]: 球: 距离圆心距离为R的点的集合. Σdim_i² = R²
*	[算法]: 计算正象限的点坐标，然后通过取负绘制其他象限.
---------------------------------------------------------------------------*/
void GraphicsND::drawSuperSphere(Mat<double>& center, double r, bool FACE, bool LINE) {
	unsigned int Dim = center.rows, maxCode = 0, times = 1, cur;
	double delta = 0.1, tmp;
	Mat<double> point(Dim, 1), tmpMat;
	for (int dim = 0; dim < Dim; dim++) { times *= 1.0 / delta + 1; maxCode += 1 << dim; }
	for (int i = 0; i < times; i++) {
		//[1] 计算正象限的点坐标
		cur = 0; point[cur] += delta;
		while (point[cur] > 1 && cur + 1 < Dim - 1) { point[cur] = 0; cur++; point[cur] += delta; }
		tmp = 1; for (int j = 0; j < Dim - 1; j++) tmp -= point[j] * point[j];
		point[Dim - 1] = sqrt(tmp);
		//[2] 然后通过取负绘制其他象限
		for (unsigned int code = 0; code <= maxCode; code++) {
			tmpMat = point;
			for (int j = 0; j < Dim; j++)
				if (code & (1 << j))
					tmpMat[j] *= -1;
			drawPoint(tmpMat.add(center, tmpMat.mult(r, tmpMat)));
		}
	}
}
/*--------------------------------[ 画网格 ]--------------------------------
*	[过程]:
		[1] 计算每一个格点的坐标
		[2] 绘制该格点对应的, 各维度方向的从min[dim] -> max[dim]的直线段
---------------------------------------------------------------------------*/
void GraphicsND::drawGrid(Mat<double>& delta, Mat<double>& max, Mat<double>& min, bool LINE) {
	int times = 1, cur;
	for (int dim = 0; dim < min.rows; dim++) times *= (max[dim] - min[dim]) / delta[dim] + 1;

	Mat<double> point(min), st, ed; point[0] -= delta[0];
	for (int i = 0; i < times; i++) {
		//[1]
		cur = 0; point[cur] += delta[cur];
		while (point[cur] > max[cur]) { point[cur] = min[cur]; cur++; point[cur] += delta[cur]; }
		//[2]
		if (!LINE)drawPoint(point);
		if (LINE) {
			st = point; ed = point;
			for (int dim = 0; dim < min.rows; dim++) {
				st[dim] = min[dim]; ed[dim] = max[dim];
				drawLine(st, ed);
				st[dim] = point[dim]; ed[dim] = point[dim];
			}
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
/*---------------- 色谱 ----------------*/
ARGB GraphicsND::colorlist(double index, int model)
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

/******************************************************************************

*                    Transformation - 3D

******************************************************************************/
/*--------------------------------[ 平移 ]--------------------------------
[1 ]   [ 1  0  0  0 ] [1]
[x'] = [dx  1  0  0 ] [x]
|y'|   |dy  0  1  0 | |y|
|z'|   |dz  0  0  1 | |z|
**-----------------------------------------------------------------------*/
void GraphicsND::translation(Mat<double>& delta, Mat<double>& transMat) {
	Mat<double> translationMat(transMat.rows);
	for (int i = 0; i < delta.rows; i++) translationMat(i + 1, 0) = delta[i];
	transMat.mult(translationMat, transMat);
}
/*--------------------------------[ 旋转 ]--------------------------------
*	[公式]: v' = q v q`¹
		q = [cos(θ/2), u sin(θ/2)]
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
	Mat<double> tmp, rotateMat(transMat.rows);
	translation(center.negative(tmp), transMat);
	// q
	if (transMat.rows - 1 == 2) { //二维 S02
		double t[] = {
			1, 0, 0,
			0, cos(theta), -sin(theta),
			0, sin(theta),  cos(theta),
		}; rotateMat.getData(t);
	}
	if (transMat.rows - 1 == 3) { //三维 S03·四元数
		Mat<double> q(transMat.rows, 1);
		rotateAxis.normalized();
		q[0] = cos(theta / 2); for (int i = 0; i < rotateAxis.rows; i++) q[i + 1] = sin(theta / 2) * rotateAxis[i];
		// rotate mat
		for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) rotateMat(i, j) = q[((j % 2 == 0 ? 1 : -1) * i + j + 4) % 4];
		for (int i = 1; i < 4; i++) rotateMat(i, i % 3 + 1) = -rotateMat(i, i % 3 + 1);
		Mat<double> rotateMatR(rotateMat);
		for (int i = 1; i < 4; i++) { rotateMat(0, i) *= -1; rotateMatR(i, 0) *= -1; }
		rotateMat.mult(rotateMat, rotateMatR);
	}
	if (transMat.rows - 1 == 4) { //四维 S04
		double c1 = cos(rotateAxis[0]), s1 = sin(rotateAxis[0]);
		double c2 = cos(rotateAxis[1]), s2 = sin(rotateAxis[1]);
		double t[] = {
			1,0,0,0,0,
			0,c1,-s1, c2,-s2,
			0,s1, c1, s1, c1,
			0,c2,-s1, c2,-s2,
			0,s2, c1, s2, c2
		}; rotateMat.getData(t);
	}
	// rotate
	transMat.mult(rotateMat, transMat);
	translation(center, transMat);
}
/*--------------------------------[ 缩放 ]--------------------------------
[1 ]   [ 1  0  0  0 ] [1]
[x'] = [ 0 sx  0  0 ] [x]
|y'|   | 0  0 sy  0 | |y|
|z'|   | 0  0  0 sz | |z|
**-----------------------------------------------------------------------*/
void GraphicsND::scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat) {
	Mat<double> tmp;
	translation(center.negative(tmp), transMat);
	// scaling
	Mat<double> scaleMat(transMat.rows);
	for (int i = 0; i < scale.rows; i++)scaleMat(i + 1, i + 1) = scale[i];
	transMat.mult(scaleMat, transMat);
	translation(center, transMat);
}