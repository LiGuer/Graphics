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
#include "RayTracing.h"
/*--------------------------------[ 初始化 ]--------------------------------*/
void RayTracing::init(int width, int height) {
	g.init(width, height);
}
/*--------------------------------[ 渲染 ]--------------------------------
*	[过程]:
		[1] 计算屏幕矢量、屏幕X,Y向轴
		[2] 对屏幕每个像素遍历
			[3] 计算像素矢量、光线矢量、光线追踪起点
			[4] 光线追踪算法
			[5] 基于结果绘制该像素色彩
-------------------------------------------------------------------------*/
void RayTracing::paint() {
	//[1]
	Mat<double> ScreenVec, ScreenXVec, ScreenYVec(3, 1);
	ScreenVec.add(gCenter, Eye.negative(ScreenVec));												//屏幕轴由眼指向屏幕中心
	{ double t[] = { ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0],1,0 }; ScreenYVec.getData(t).normalization(); }	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalization();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	//[2]
	double minDistance = 0, RayFaceDistance;
	Mat<double> PixYVec, PixXVec, PixVec, Ray, RaySt;
	for (int x = 0; x < g.Canvas.rows; x++) {
		for (int y = 0; y < g.Canvas.cols; y++) {
			//[3]
			PixVec.add(PixXVec.mult(x - g.Canvas.rows / 2, ScreenXVec), PixYVec.mult(y - g.Canvas.cols / 2, ScreenYVec));
			Ray.add(ScreenVec, PixVec).normalization();
			RaySt.add(gCenter, PixVec);
			//[4][5]
			RGB color(0);
			traceRay(RaySt, Ray, color, 0);
			g.setPoint(g.Canvas.rows - x, y, color.R * 0x10000 + color.G * 0x100 + color.B);
		}
	}
}
/*--------------------------------[ 追踪光线 ]--------------------------------
*	[过程]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
-----------------------------------------------------------------------------*/
RayTracing::RGB RayTracing::traceRay(Mat<double>& RaySt, Mat<double>& Ray, RGB& color, int level) {
	double minDistance = DBL_MAX;
	Mat<double> intersection, intersectionTmp, FaceVec, FaceVecTmp;
	Triangle closestTriangle;
	//[1]
	for (int i = 0; i < TriangleSet.size(); i++) {
		//[2][3]
		double distance = seekIntersection(TriangleSet[i], RaySt, Ray, FaceVecTmp, intersectionTmp);
		if (distance > 0 && distance < minDistance) {
			minDistance = distance; closestTriangle = TriangleSet[i];
			FaceVec = FaceVecTmp; intersection = intersectionTmp;
		}
	}
	//[4]
	if (minDistance != DBL_MAX && level < maxRayLevel) {
		// Reflex Ray
		double RayFaceCosAngle = FaceVec.dot(Ray);
		if (RayFaceCosAngle > 0) FaceVec.negative(FaceVec);
		Mat<double> Reflect;
		Reflect.add(Reflect.mult(pow(1 / RayFaceCosAngle, 2), FaceVec), Ray).normalization();
		traceRay(intersection, Reflect, color, level + 1);
		// Reflex Rate
		double k = closestTriangle.material->reflexRate;
		color.R *= k * closestTriangle.material->color.R / 255.0;
		color.G *= k * closestTriangle.material->color.G / 255.0;
		color.B *= k * closestTriangle.material->color.B / 255.0;
	}
	else {
		Mat<double> Light;
		for (int i = 0; i < LightSource.size(); i++) {
			Light.add(LightSource[i], Ray.negative(Light));
			double LightSourceAngle = (Ray.dot(Light) / Light.norm() + 1) / 2;
			color.R = 0xFF * LightSourceAngle;
			color.G = 0xFF * LightSourceAngle;
			color.B = 0xFF * LightSourceAngle;
		}
	}
	return color;
}
/*--------------------------------[ 求交点 ]--------------------------------
*	[流程]:
		[1] 计算三角形所在面矢量
		[2] 计算光线面交点、光线面相交所走过距离
		[3] 判断交点是否在三角形内部, 若否返回-1
*	[算法]:
		平面方程: Af (X - Xf) + BY (Y - Yf) + C (Z - Zf) = 0
		直线方程: (X - Xl) / Al = (Y - Yl) / Bl = (Z - Zl) / Cl = K
		点面距:	d = |AXp + BYp + CZp + D| / sqrt(A² + B² + C²)
		线面交点: K = [Af(Xf - Xl) + Bf(Yf - Yl) + Cf(Zf - Zl)] / (Af Al + Bf Bl + Cf Cl) 即光线走过线距离
				  X = K Al + Xl
---------------------------------------------------------------------------*/
double RayTracing::seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, Mat<double>& intersection) {
	//[1]
	Mat<double> edge[2], tmp;
	edge[0].add(triangle.p[1], triangle.p[0].negative(edge[0]));
	edge[1].add(triangle.p[2], triangle.p[0].negative(edge[1]));
	FaceVec.crossProduct(edge[0], edge[1]).normalization();
	//[2]
	double RayFaceDistance = FaceVec.dot(tmp.add(triangle.p[0], RaySt.negative(tmp))) / FaceVec.dot(Ray);
	intersection.add(tmp.mult(RayFaceDistance, Ray), RaySt);
	//[3]
	Mat<double> tmpEdge; tmpEdge.add(intersection, triangle.p[0].negative(tmpEdge));
	double inverDeno = 1 / (edge[0].dot(edge[0]) * edge[1].dot(edge[1]) - edge[0].dot(edge[1]) * edge[1].dot(edge[0]));
	double u = (edge[1].dot(edge[1]) * edge[0].dot(tmpEdge) - edge[0].dot(edge[1]) * edge[1].dot(tmpEdge)) * inverDeno;
	double v = (edge[0].dot(edge[0]) * edge[1].dot(tmpEdge) - edge[1].dot(edge[0]) * edge[0].dot(tmpEdge)) * inverDeno;						
	if (u < 0 || u > 1 || v < 0 || v > 1 || u + v > 1) return -1;		// if u,v out of range, return directly
	return  RayFaceDistance;
}
/*--------------------------------[ 画三角形 ]--------------------------------*/
void RayTracing::drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Material* material) {
	Triangle triangle;
	triangle.p[0] = p1;
	triangle.p[1] = p2;
	triangle.p[2] = p3;
	triangle.material = material;
	TriangleSet.push_back(triangle);
}
/*--------------------------------[ 画矩形 ]--------------------------------*/
void RayTracing::drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct) {

}
/*--------------------------------[ 画四边形 ]--------------------------------*/
void RayTracing::drawQuadrilateral(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, Material* material) {
		drawTriangle(p1, p2, p3, material); 
		drawTriangle(p1, p4, p3, material);
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
void RayTracing::drawPolygon(Mat<double> p[], int n, Material* material) {
	for (int k = 1; k <= (n + 2) / 3; k++)
		for (int i = 0; i <= n - 2 * k; i += 2 * k)
			drawTriangle(p[i], p[i + k], p[(i + 2 * k) % n], material);
}
/*--------------------------------[ 画圆 ]--------------------------------
**-----------------------------------------------------------------------*/
void RayTracing::drawCircle(Mat<double>& center, double r, Mat<double>* direct) {

}
/*--------------------------------[ 画椭圆 ]--------------------------------
**-----------------------------------------------------------------------*/
void RayTracing::drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct) {
}
/*--------------------------------[ 画曲面 ]--------------------------------*/
void RayTracing::drawSurface(Mat<double> z, double xs, double xe, double ys, double ye) {
	Mat<double> p(3, 1), pl(3, 1), pu(3, 1);
	double dx = (xe - xs) / z.rows, dy = (ye - ys) / z.cols;
	for (int y = 0; y < z.cols; y++) {
		for (int x = 0; x < z.rows; x++) {
			{double t[] = { xs + x * dx,ys + y * dy,z(x,y) }; p.getData(t); }
			if (x > 0) {
				double t[] = { xs + (x - 1) * dx,ys + y * dy,z(x - 1,y) };
				//pl.getData(t); drawLine(pl, p);
			}
			if (y > 0) {
				double t[] = { xs + x * dx,ys + (y - 1) * dy,z(x,y - 1) };
				//pu.getData(t); drawLine(pu, p);
			}
		}
	}
}
/*--------------------------------[ 画四面体 ]--------------------------------*/
void RayTracing::drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, Material* material) {
	drawTriangle(p1, p2, p3, material);
	drawTriangle(p2, p3, p4, material);
	drawTriangle(p3, p4, p1, material);
	drawTriangle(p4, p1, p2, material);
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
void RayTracing::drawCuboid(Mat<double>& pMin, Mat<double>& pMax, Material* material) {
	Mat<double> pMinTmp[3], pMaxTmp[3];
	for (int i = 0; i < 3; i++) {
		pMinTmp[i] = pMin; pMinTmp[i][i] = pMax[i];
		pMaxTmp[i] = pMax; pMaxTmp[i][i] = pMin[i];
	}
	for (int i = 0; i < 3; i++) {
		drawTriangle(pMin, pMinTmp[i], pMaxTmp[(i + 2) % 3], material);
		drawTriangle(pMax, pMaxTmp[i], pMinTmp[(i + 2) % 3], material);
		drawTriangle(pMin, pMinTmp[(i + 1) % 3], pMaxTmp[(i + 2) % 3], material);
		drawTriangle(pMax, pMaxTmp[(i + 1) % 3], pMinTmp[(i + 2) % 3], material);
	}
}
/*--------------------------------[ 画圆台 ]--------------------------------
* [过程]:
		[1] 计算旋转矩阵
		[2] 根据旋转矩阵, 计算绘制点坐标, 完成绘制
**------------------------------------------------------------------------*/
void RayTracing::drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta, Material* material) {
	// 计算 Rotate Matrix
	Mat<double> direction, rotateAxis, rotateMat(4), zAxis(3, 1), tmp; {double t[] = { 0, 0, 1 }; zAxis.getData(t); }
	direction.add(ed, st.negative(direction));
	if (direction[0] != 0 || direction[1] != 0) {
		//rotate(
		//	rotateAxis.crossProduct(direction, zAxis),
		//	-acos(tmp.dot(direction, zAxis) / direction.norm()),
		//	tmp.zero(3, 1), rotateMat
		//);
		rotateMat.cut(1, 3, 1, 3, rotateMat);
	}
	else rotateMat.E(3);
	// 画圆台
	Mat<double> stPoint, edPoint, preStPoint, preEdPoint, deltaVector(3, 1);
	for (int i = 0; i <= 360 / delta; i++) {
		{double t[] = { cos(i * delta * 2.0 * PI / 360), sin(i * delta * 2.0 * PI / 360),0 }; deltaVector.getData(t); }
		deltaVector.mult(rotateMat, deltaVector);
		stPoint.add(st, stPoint.mult(Rst, deltaVector));
		edPoint.add(ed, edPoint.mult(Red, deltaVector));
		if (i != 0) {
			drawTriangle(stPoint, preStPoint, edPoint);
			drawTriangle(preStPoint, preEdPoint, edPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}
/*--------------------------------[ 画圆柱 ]--------------------------------*/
void RayTracing::drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta, Material* material) {
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
void RayTracing::drawSphere(Mat<double>& center, double r, int delta, Material* material) {
	// 经纬度法
	Mat<double> point(3, 1);
	for (int i = 0; i < 360 / delta; i++) {
		double theta = (i * delta) * 2.0 * PI / 360;
		for (int j = -90 / delta; j <= 90 / delta; j++) {
			double phi = (j * delta) * 2.0 * PI / 360;
			point[0] = r * cos(phi) * cos(theta) + center[0];
			point[1] = r * cos(phi) * sin(theta) + center[1];
			point[2] = r * sin(phi) + center[2];
		}
	}
}
/*--------------------------------[ getSphereFibonacciPoint 球面均匀点分布 ]--------------------------------
*	[Referance]:
		[1] Thanks and copyright for https://github.com/SebLague/Boids
**---------------------------------------------------------------------------------------------------------*/
void RayTracing::drawSphere2(Mat<double>& center, double r, int n, Material* material) {
	// 均匀球面点
	Mat<double> point(3, 1);
	double goldenRatio = (1 + sqrt(5)) / 2;				// 黄金分割点
	double angleIncrement = PI * 2 * goldenRatio;
	for (int i = 0; i < 300; i++) {
		double t = (double)i / n, inclination = acos(1 - 2 * t), azimuth = angleIncrement * i;
		point[0] = center[0] + r * sin(inclination) * cos(azimuth);
		point[1] = center[1] + r * sin(inclination) * sin(azimuth);
		point[2] = center[2] + r * cos(inclination);
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
void RayTracing::drawEllipsoid(Mat<double>& center, Mat<double>& r, Material* material) {
	const int delta = 5;
	Mat<double> point(3, 1);
	for (int i = 0; i < 360 / delta; i++) {
		double theta = (i * delta) * 2.0 * PI / 360;
		for (int j = -90 / delta; j <= 90 / delta; j++) {
			double phi = (j * delta) * 2.0 * PI / 360;
			point[0] = r[0] * cos(phi) * cos(theta) + center[0];
			point[1] = r[1] * cos(phi) * sin(theta) + center[1];
			point[2] = r[2] * sin(phi) + center[2];
		}
	}
}