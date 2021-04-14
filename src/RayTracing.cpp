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
	{ double t[] = { ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0],1,0 }; ScreenYVec.getData(t).normalized(); }	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalized();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	//[2]
	double minDistance = 0, RayFaceDistance;
	Mat<double> PixYVec, PixXVec, PixVec, Ray, RaySt;
	for (int x = 0; x < g.Canvas.rows; x++) {
		printf("%d\t", x);
		for (int y = 0; y < g.Canvas.cols; y++) {
			//[3]
			PixVec.add(PixXVec.mult(x - g.Canvas.rows / 2, ScreenXVec), PixYVec.mult(y - g.Canvas.cols / 2, ScreenYVec));
			Ray.add(ScreenVec, PixVec).normalized();
			RaySt.add(gCenter, PixVec);
			//[4][5]
			RGB color(0);
			traceRay(RaySt, Ray, color, 0);
			g.setPoint(g.Canvas.rows - x, y, color.R * 0x10000 + color.G * 0x100 + color.B);
		}
	}
}
/*--------------------------------[ 追踪光线 ]--------------------------------
*	[算法]:
			设面矢F, 入射光L
			[反射光]: 反射定律: 入射角 == 反射角
					  Lf = L - F·2 cos<L,F>
			[折射光]: 折射定律: n1·sinθ1 = n2·sinθ2 
					  => sin α =  n·sin <L,F>  ,  n = n入 / n折
					  Lz = L + F·sin(α - <L,F>) / sinα
					     = L + F·( cos<L,F> - sqrt( 1/n² - 1 + cos²<L,F> ) )
*	[过程]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
-----------------------------------------------------------------------------*/
RGB RayTracing::traceRay(Mat<double>& RaySt, Mat<double>& Ray, RGB& color, int level) {
	double minDistance = DBL_MAX;
	Mat<double> intersection, intersectionTmp, FaceVec, FaceVecTmp;
	Material* intersectMaterial = NULL;
	//[1]
	for (int i = 0; i < TriangleSet.size(); i++) {
		//[2][3]
		double distance = seekIntersection(TriangleSet[i], RaySt, Ray, FaceVecTmp, intersectionTmp);
		if (distance > 0.1 && distance < minDistance) {		//distance > 1而不是> 0，是因为反射光线在接触面的精度内，来回碰自己....
			minDistance = distance; intersectMaterial = TriangleSet[i].material;
			FaceVec = FaceVecTmp; intersection = intersectionTmp;
		}
	}
	//[4]
	if (minDistance != DBL_MAX && level < maxRayLevel && intersectMaterial != NULL) {
		Mat<double> RayTmp;
		RGB colorTmp(0);
		if (intersectMaterial->rediateRate != 0) {					//Light Source
			return color = intersectMaterial->color;
		}
		if (intersectMaterial->diffuseReflect != 0) {				//diffuseReflect
			double LightSourceAngle = 0;
			for (int i = 0; i < LightSource.size(); i++) {
				RayTmp.add(LightSource[i], intersection.negative(RayTmp)).normalized();
				double LightSourceAngleTmp = (FaceVec.dot(Ray) > 0 ? -1 : 1)* FaceVec.dot(RayTmp);
				LightSourceAngle = LightSourceAngle > LightSourceAngleTmp ? LightSourceAngle : LightSourceAngleTmp;
			}
			colorTmp = 0xFFFFFF;
			colorTmp *= LightSourceAngle;
			goto TRACERAY_LABLE;
			//ColorBlend_Add(color, colorTmp, color);
		}
		if (intersectMaterial->reflectRate != 0) {				//reflect
			traceRay(intersection, reflect(Ray, FaceVec, RayTmp), colorTmp = 0, level + 1);
			colorTmp *= intersectMaterial->reflectRate;
			ColorBlend_Add(color, colorTmp, color);
		}
		if (intersectMaterial->refractRate != 0) {				//refract
			double refractRateBufferTmp = refractRateBuffer; refractRateBuffer = refractRateBuffer == intersectMaterial->refractRate ? 1 : intersectMaterial->refractRate;
			refract(Ray, FaceVec, RayTmp, refractRateBufferTmp, refractRateBuffer);
			intersectionTmp.add(intersectionTmp.mult(0.1, RayTmp), intersection);
			traceRay(intersectionTmp, RayTmp, colorTmp = 0, level + 1);
			refractRateBuffer = refractRateBufferTmp;
			colorTmp *= 1 - intersectMaterial->reflectRate;		//###菲涅耳方程 未完成
			ColorBlend_Add(color, colorTmp, color);
			
		}
		TRACERAY_LABLE:
		color.R *= (double)intersectMaterial->color.R / 0xFF;
		color.G *= (double)intersectMaterial->color.G / 0xFF;
		color.B *= (double)intersectMaterial->color.B / 0xFF;
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
		球方程: (X - Xs)² + (Y - Ys)² + (Z - Zs)² = R²
		线球交点: K²(Al² + Bl² + Cl²) + 2K(Al ΔX + Bl ΔY + Cl ΔZ) + (ΔX² + ΔY² + ΔZ² - R²) = 0
				  ΔX = Xl - Xs
				  Δ = b² - 4ac = 4(Al ΔX + Bl ΔY + Cl ΔZ)² - 4(Al² + Bl² + Cl²)(ΔX² + ΔY² + ΔZ² - R²)
				  若Δ≥0 有交点.
				  K = ( -b ± sqrt(Δ) ) / 2a	即光线走过线距离
---------------------------------------------------------------------------*/
double RayTracing::seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, Mat<double>& intersection) {
	// Sphere Seek Intersection
	if (triangle.p[2][0] == NULL) {
		// 计算ΔX、Δ
		double R = triangle.p[1][0], Delta;
		Mat<double> RayStCenter;
		RayStCenter.add(RaySt, triangle.p[0].negative(RayStCenter));
		Delta = 4 * pow(Ray.dot(RayStCenter), 2) - 4 * Ray.dot(Ray) * (RayStCenter.dot(RayStCenter) - R * R);
		if (Delta < 0) return -1;
		// 有交点，计算RayFaceDistance、intersection、FaceVec
		double RayFaceDistance = -2 * Ray.dot(RayStCenter);
		RayFaceDistance += RayFaceDistance - sqrt(Delta) > 0 ? -sqrt(Delta) : +sqrt(Delta);
		RayFaceDistance /= 2 * Ray.dot(Ray);
		intersection.add(intersection.mult(RayFaceDistance, Ray), RaySt);
		FaceVec.add(intersection, triangle.p[0].negative(FaceVec)).normalized();
		return RayFaceDistance;
	}
	//[1]
	Mat<double> edge[2], tmp;
	edge[0].add(triangle.p[1], triangle.p[0].negative(edge[0]));
	edge[1].add(triangle.p[2], triangle.p[0].negative(edge[1]));
	FaceVec.crossProduct(edge[0], edge[1]).normalized();
	//[2]
	double RayFaceDistance = FaceVec.dot(tmp.add(triangle.p[0], RaySt.negative(tmp))) / FaceVec.dot(Ray);
	intersection.add(intersection.mult(RayFaceDistance, Ray), RaySt);
	//[3]
	Mat<double> tmpEdge; tmpEdge.add(intersection, triangle.p[0].negative(tmpEdge));
	double inverDeno = 1 / (edge[0].dot(edge[0]) * edge[1].dot(edge[1]) - edge[0].dot(edge[1]) * edge[1].dot(edge[0]));
	double u = (edge[1].dot(edge[1]) * edge[0].dot(tmpEdge) - edge[0].dot(edge[1]) * edge[1].dot(tmpEdge)) * inverDeno;
	double v = (edge[0].dot(edge[0]) * edge[1].dot(tmpEdge) - edge[1].dot(edge[0]) * edge[0].dot(tmpEdge)) * inverDeno;						
	if (u < 0 || u > 1 || v < 0 || v > 1 || u + v > 1) return -1;		// if u,v out of range, return directly
	return RayFaceDistance;
}
/*--------------------------------[ 反射 ]--------------------------------*/
Mat<double>& RayTracing::reflect(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& reflectRay) {
	return reflectRay.add(reflectRay.mult(-2 * faceVec.dot(incidentRay), faceVec), incidentRay).normalized();
}
/*--------------------------------[ 折射 ]--------------------------------*/
Mat<double>& RayTracing::refract(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& refractRay, double rateIn, double rateOut) {
	double refractRate = rateIn / rateOut;
	double CosIn = faceVec.dot(incidentRay);
	double CosOut = 1 - refractRate * refractRate * (1 - CosIn * CosIn);
	if (CosOut >= 0) reflect(incidentRay, faceVec, refractRay);			//全反射
	Mat<double> tmp;
	CosOut = sqrt(CosOut);
	return refractRay.add(refractRay.mult(refractRate, incidentRay), tmp.mult((CosIn > 0 ? 1 : -1)* CosOut - refractRate * CosIn, faceVec)).normalized();
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
		GraphicsND::rotate(
			rotateAxis.crossProduct(direction, zAxis),
			-acos(tmp.dot(direction, zAxis) / direction.norm()),
			tmp.zero(3, 1), rotateMat
		);
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
			drawTriangle(stPoint, preStPoint, edPoint, material);
			drawTriangle(preStPoint, preEdPoint, edPoint, material);
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
**-----------------------------------------------------------------------*/
void RayTracing::drawSphere(Mat<double>& center, double r, Material* material) {
	Mat<double> p1(3, 1), p2 = (3, 1);
	p1[0] = p1[1] = p1[2] = r; p2[0] = NULL;
	drawTriangle(center, p1, p2, material); 
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