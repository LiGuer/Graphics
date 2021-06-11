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
[Reference]:
[1] Thanks for Kevin Beason at http://www.kevinbeason.com/smallpt/
==============================================================================*/
#include "RayTracing.h"
#include "GraphicsFileCode.h"
#include <time.h>
/*--------------------------------[ 初始化 ]--------------------------------*/
void RayTracing::init(int width, int height) {
	ScreenPix.zero(height, width);
	Screen.   zero(height, width);
	for (int i = 0; i < Screen.size(); i++) Screen[i].zero(3);
}
/*--------------------------------[ 画像素 ]--------------------------------*/
void RayTracing::setPix(int x, int y, Mat<>& color) {
	if (x < 0 || x >= ScreenPix.rows || y < 0 || y >= ScreenPix.cols) return;
	ScreenPix(ScreenPix.rows - x - 1, y).R = std::min((int)(color[0] * 0xFF), 0xFF);
	ScreenPix(ScreenPix.rows - x - 1, y).G = std::min((int)(color[1] * 0xFF), 0xFF);
	ScreenPix(ScreenPix.rows - x - 1, y).B = std::min((int)(color[2] * 0xFF), 0xFF);
}
/*--------------------------------[ 渲染 ]--------------------------------
*	[过程]:
		[1] 计算屏幕矢量、屏幕X,Y向轴
		[2] 对屏幕每个像素遍历
			[3] 计算像素矢量、光线矢量、光线追踪起点
			[4] 光线追踪算法
			[5] 基于结果绘制该像素色彩
-------------------------------------------------------------------------*/
void RayTracing::paint(const char* fileName, int sampleSt, int sampleEd) {
	//[1]
	static Mat<> ScreenVec, ScreenXVec, ScreenYVec(3);
	ScreenVec. sub(gCenter, Eye);																	//屏幕轴由眼指向屏幕中心
	ScreenYVec.getData(ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0], 1, 0).normalized();	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalized();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	//[2]
	static Mat<> PixYVec, PixXVec, PixVec, Ray, RaySt, color(3);
	for (int sample = sampleSt; sample < sampleEd; sample++) {
		GraphicsFileCode::ppmWrite(fileName, ScreenPix); clock_t start = clock();
		double rate = 1.0 / (sample + 1);
		for (int x = 0; x < Screen.rows; x++) {
			for (int y = 0; y < Screen.cols; y++) {
				PixVec.add(														//[3]
					PixXVec.mul(x + RAND_DBL - Screen.rows / 2 - 0.5, ScreenXVec),
					PixYVec.mul(y + RAND_DBL - Screen.cols / 2 - 0.5, ScreenYVec)
				);
				traceRay(														//[4][5]
					RaySt.add(gCenter,   PixVec), 
					Ray.  add(ScreenVec, PixVec).normalized(), 
					color.zero(), 0
				);
				setPix(x, y, (Screen(x, y) *= 1 - rate) += (color *= rate));
			} 
		} printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));
	}
}
/******************************************************************************
*						追踪光线
*	[算法]:
*	[过程]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
&	[注]:distance > 1而不是> 0，是因为反射光线在接触面的精度内，来回碰自己....
******************************************************************************/
Mat<>& RayTracing::traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level) {
	//[1][2][3]
	double minDis = DBL_MAX, dis;
	Triangle* minDisTri = NULL;
	for (int i = 0; i < TriangleSet.size(); i++) {
		dis    = seekIntersection(TriangleSet[i], RaySt, Ray);
		minDis = (dis > eps && dis < minDis) ? minDisTri = &TriangleSet[i], dis : minDis;
	}
	if (minDis == DBL_MAX)			return color.zero();			//Miss intersect
	Material* material = minDisTri->material;
	if (material->rediateRate != 0)	return color = material->color;	//Light Source
	if (level > maxRayLevel)		return color.zero();			//Max Ray Level
	//[4] intersection & FaceVec
	Mat<> intersection, faceVec, tmp, rayTmp;
	{
		intersection.add(RaySt, intersection.mul(minDis, Ray));
		if (minDisTri->p[2][0] == DBL_MAX)
			 faceVec.sub(intersection, minDisTri->p[0]).normalized();
		else faceVec.crossProduct_(
				   tmp.sub(minDisTri->p[1], minDisTri->p[0]), 
				rayTmp.sub(minDisTri->p[2], minDisTri->p[0])
			 ).normalized();
	}
	//[5]
	if (material->quickReflect != 0) {						//Quick Reflect: 若该点处的表面是(快速)散射面，计算点光源直接照射该点产生的颜色
		double lightCos = 0, t;
		faceVec *= faceVec.dot(Ray) > 0 ? -1 : 1;
		for (int i = 0; i < PointLight.size(); i++) {
			t = (faceVec.dot(rayTmp.sub(PointLight[i], intersection).normalized()) + 1) / 2;
			lightCos = t > lightCos ? t : lightCos;
		} color.ones(3) *= lightCos;
	}
	else if (material->diffuseReflect != 0) {				//Diffuse Reflect
		traceRay(intersection, diffuseReflect(Ray, faceVec, rayTmp),	color.zero(), level + 1);
		color *= material->reflectRate;
	}
	else if (material->refractRate != 0) {					//Refract
		double t = refractRateBuf; refractRateBuf = refractRateBuf == material->refractRate ? 1 : material->refractRate;
		refract(Ray, faceVec, rayTmp, t, refractRateBuf);
		traceRay(tmp.add(tmp.mul(eps, rayTmp), intersection), rayTmp,	color.zero(), level + 1);
		refractRateBuf = t;
	}
	else if (material->reflectRate != 0) {					//Reflect
		traceRay(intersection, reflect(Ray, faceVec, rayTmp),			color.zero(), level + 1);
		color *= material->reflectRate;
	}
	return color.elementMul(material->color);
}
/******************************************************************************
						求交点
*	[算法]:
		[射线球面交点]
		球方程: (X - Xs)² + (Y - Ys)² + (Z - Zs)² = R²
		线球交点: K²(Al² + Bl² + Cl²) + 2K(Al ΔX + Bl ΔY + Cl ΔZ) + (ΔX² + ΔY² + ΔZ² - R²) = 0
				  ΔX = Xl - Xs
				  Δ = b² - 4ac = 4(Al ΔX + Bl ΔY + Cl ΔZ)² - 4(Al² + Bl² + Cl²)(ΔX² + ΔY² + ΔZ² - R²)
				  若Δ≥0 有交点.
				  K = ( -b ± sqrt(Δ) ) / 2a	即光线走过线距离
		[射线三角形交点]
		Moller-Trumbore方法(1997)
		射线: P = O + t D
		射线三角形交点:
				O + t D = (1 - u - v)V0 + u V1 + v V2
				[ -D  V1-V0  V2-V0] [ t  u  v ]' = O - V0
				T = O - V0    E1 = V1 - V0    E2 = V2 - V0
				[ -D  E1  E2 ] [ t  u  v ]' = T
				t = | T  E1  E2| / |-D  E1  E2|
				u = |-D   T  E2| / |-D  E1  E2|
				v = |-D  E1  E2| / |-D  E1  E2|
		(混合积公式): |a  b  c| = a×b·c = -a×c·b
				t = (T×E1·E2) / (D×E2·E1)
				u = (D×E2· T) / (D×E2·E1)
				v = (T×E1· D) / (D×E2·E1)
******************************************************************************/
double RayTracing::seekIntersection(Triangle& triangle, Mat<>& RaySt, Mat<>& Ray) {
	if (triangle.p[2][0] == DBL_MAX) { // Sphere Seek Intersection
		static Mat<> RayStCenter; RayStCenter.sub(RaySt, triangle.p[0]);
		double R = triangle.p[1][0], 
			   A = Ray.dot(Ray), 
			   B = 2 * Ray.dot(RayStCenter),
		       Delta = B * B - 4 * A * (RayStCenter.dot(RayStCenter) - R * R);
		if (Delta < 0) return -DBL_MAX;									//有无交点
		Delta = sqrt(Delta);
		return (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
	}
	else { // Triangle Seek Intersection
		static Mat<> edge[2], tmp, p, q;
		edge[0].sub(triangle.p[1], triangle.p[0]);
		edge[1].sub(triangle.p[2], triangle.p[0]);
		// p & a & tmp
		static double a, u, v;
		a = p.crossProduct_(Ray, edge[1]).dot(edge[0]);
		if (a > 0) tmp.sub(RaySt, triangle.p[0]);
		else       tmp.sub(triangle.p[0], RaySt), a = -a;
		if (a < 1e-4)		return -DBL_MAX;								//射线与三角面平行
		// u & q & v
		u = p.dot(tmp) / a;
		if(u < 0 || u > 1)	return -DBL_MAX;
		v = q.crossProduct_(tmp, edge[0]).dot(Ray) / a;
		return (v < 0 || u + v > 1) ? -DBL_MAX : q.dot(edge[1]) / a;
	}
}
/*--------------------------------[ 画三角形 ]--------------------------------*/
void RayTracing::drawTriangle(Mat<>& p1, Mat<>& p2, Mat<>& p3, Material* material) {
	Triangle triangle;
	triangle.p[0] = p1;
	triangle.p[1] = p2;
	triangle.p[2] = p3;
	triangle.material = material;
	TriangleSet.push_back(triangle); 
}
/*--------------------------------[ 画球 ]--------------------------------*/
void RayTracing::drawSphere(Mat<>& center, double r, Material* material) {
	Mat<> p1(3), p2(3);
	p1[0] = p1[1] = p1[2] = r; p2[0] = DBL_MAX;
	drawTriangle(center, p1, p2, material); 
}
/*#############################################################################

*						几何光学  Geometrical Optics

##############################################################################*/
/******************************************************************************
*						反射
*	[反射定律]: 入射角 == 反射角
			设面矢F, 入射光L
			Lf = L - F·2 cos<L,F>
******************************************************************************/
Mat<>& RayTracing::reflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
	return RayO.add(RayO.mul(-2 * faceVec.dot(RayI), faceVec), RayI).normalized();
}
/******************************************************************************
*						折射
*	[折射定律]: n1·sinθ1 = n2·sinθ2
			设面矢F, 入射光L
			=> sin θo =  n·sin θi  ,  n = n入 / n折, θi = <L,F>
			Lo = Li + F·sin(θo - θi) / sinθo
			   = Li + F·( cosθi·sinθo - cosθo·sinθi) / sinθo
			   = Li + F·( cosθi - cosθo / n )
			   = Li + F·( cosθi - sqrt(1 - n²·sin²θi) / n )
******************************************************************************/
Mat<>& RayTracing::refract(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO) {
	double k    = rateI / rateO,
		   CosI = faceVec.dot(RayI),
		   CosO = 1 - pow(k, 2) * (1 - pow(CosI, 2));
	if (CosO < 0) return reflect(RayI, faceVec, RayO);			//全反射
	return RayO.add(RayI, RayO.mul(-CosI - (CosI > 0 ? -1 : 1) * sqrt(CosO) / k, faceVec)).normalized();
}
/******************************************************************************
*						漫反射
*	[算法]: MonteCarlo: 随机持续采样
*	[漫反射]: 在面矢半球内，面积均匀的随机取一射线，作为反射光线.
******************************************************************************/
Mat<>& RayTracing::diffuseReflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
	double r1 = 2 * PI * RAND_DBL, r2 = RAND_DBL;
	static Mat<> tmp(3), tmp1(3), tmp2(3);
	faceVec *= faceVec.dot(RayI) > 0 ? -1 : 1;
	tmp[0] = fabs(faceVec[0]) > 0.1 ? 0 : 1;
	tmp[1] = tmp[0] == 0 ? 1 : 0;
	tmp1.mul(cos(r1) * sqrt(r2), tmp1.crossProduct_(tmp, faceVec).normalized());
	tmp2.mul(sin(r1) * sqrt(r2), tmp2.crossProduct_(faceVec, tmp1));
	return RayO.add(RayO.add(RayO.mul(sqrt(1 - r2), faceVec), tmp1), tmp2).normalized();
}
