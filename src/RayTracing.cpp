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
	{ double t[] = { 1,ScreenVec[1] == 0 ? 0 : -ScreenVec[0] / ScreenVec[1],0 }; ScreenYVec.getData(t).normalization(); }	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalization();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	//[2]
	double minDistance = 0, RayFaceDistance;
	Mat<double> PixYVec, PixXVec, PixVec, Ray, RaySt;
	for (int x = 0; x < g.Canvas.rows; x++) {
		for (int y = 0; y < g.Canvas.cols; y++) {
			//[3]
			PixVec.add(PixXVec.mult(x - g.Canvas.rows / 2, ScreenXVec), PixYVec.mult(y - g.Canvas.cols / 2, ScreenYVec));
			Ray.add(ScreenVec, PixVec);
			RaySt.add(gCenter, PixVec);
			//[4][5]
			unsigned int color = traceRay(RaySt, Ray, 0);
			g.setPoint(x, y, color);
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
unsigned int RayTracing::traceRay(Mat<double>& RaySt, Mat<double>& Ray, int level) {
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
	unsigned int color = 0;
	if (closestTriangle.material != NULL && closestTriangle.material->color != 0)
		return closestTriangle.material->color;
	if (minDistance != DBL_MAX && level < maxRayLevel) {
		Mat<double> Reflect;
		Reflect.add(Ray, Reflect.add(FaceVec,Ray.negative(Reflect)));
		color = traceRay(intersection, Reflect, level + 1);
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
		线面交点: K = [Af(Xf - Xl) + Bf(Yf - Yl) + Cf(Zf - Zl)] / (Af Al + Bf Bl + Cf Cl)
				  X = K Al + Xl
---------------------------------------------------------------------------*/
double RayTracing::seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, Mat<double>& intersection) {
	//[1]
	Mat<double> edge[2], tmp;
	edge[0].add(triangle.p[1], triangle.p[0].negative(edge[0]));
	edge[1].add(triangle.p[2], triangle.p[0].negative(edge[1]));
	FaceVec.crossProduct(edge[0], edge[1]);
	//[2]
	double K = FaceVec.dot(tmp.add(triangle.p[0], RaySt.negative(tmp))) / FaceVec.dot(Ray);
	intersection.add(tmp.mult(K, Ray), RaySt);
	double RayFaceDistance = tmp.norm();
	//[3]
	Mat<double> tmpEdge; tmpEdge.add(intersection, triangle.p[0].negative(tmpEdge));
	double inverDeno = 1 / (edge[0].dot(edge[0]) * edge[1].dot(edge[1]) - edge[0].dot(edge[1]) * edge[1].dot(edge[0]));
	double u = (edge[1].dot(edge[1]) * edge[0].dot(tmpEdge) - edge[0].dot(edge[1]) * edge[1].dot(tmpEdge)) * inverDeno;
	double v = (edge[0].dot(edge[0]) * edge[1].dot(tmpEdge) - edge[1].dot(edge[0]) * edge[0].dot(tmpEdge)) * inverDeno;						
	if (u < 0 || u > 1 || v < 0 || v > 1 || u + v > 1) return -1;		// if u,v out of range, return directly
	return  RayFaceDistance;
}