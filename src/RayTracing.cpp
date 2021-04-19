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
#include <time.h>
/*--------------------------------[ 初始化 ]--------------------------------*/
void RayTracing::init(int width, int height) {
	ScreenPix.zero(height, width);
	Screen.zero(height, width);
	for (int i = 0; i < Screen.size(); i++)Screen[i].zero(3, 1);
}
/*--------------------------------[ 读图 ]--------------------------------*/
void RayTracing::readImg(const char* fileName) {
	FILE* fin = fopen(fileName, "rb");
	int rows, cols;
	fscanf(fin, "P6\n%d %d\n255\n", &cols, &rows);						// 读图片格式、宽高、最大像素值
	init(cols, rows);
	fread(&ScreenPix[0], 1, ScreenPix.size() * 3, fin);					// 读RGB数据
	for (int x = 0; x < Screen.rows; x++)
		for (int y = 0; y < Screen.cols; y++)
			for (int k = 0; k < 3; k++)
				Screen(x, y)[k] = (double)ScreenPix(Screen.rows - 1 - x, y)[k] / 0xFF;
	fclose(fin);
}
/*--------------------------------[ 存图 ]--------------------------------*/
void RayTracing::writeImg(const char* filename) {
	FILE* fp = fopen(filename, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", ScreenPix.cols, ScreenPix.rows);	// 写图片格式、宽高、最大像素值
	fwrite(ScreenPix.data, 1, ScreenPix.size() * 3, fp);				// 写RGB数据
	fclose(fp);
}
/*--------------------------------[ 画像素 ]--------------------------------*/
void RayTracing::setPix(int x, int y, Mat<double>& color) {
	if (x < 0 || x >= ScreenPix.rows || y < 0 || y >= ScreenPix.cols) return;
	ScreenPix(ScreenPix.rows - x - 1, y).R = min((int)(color[0] * 0xFF), 0xFF);
	ScreenPix(ScreenPix.rows - x - 1, y).G = min((int)(color[1] * 0xFF), 0xFF);
	ScreenPix(ScreenPix.rows - x - 1, y).B = min((int)(color[2] * 0xFF), 0xFF);
}
/*--------------------------------[ 读3D模型 ]--------------------------------*/
void RayTracing::readObj(const char* fileName,Mat<double>& origin, Material* material) {
	FILE* fin = fopen(fileName, "rb");
	char Index = 0; double R;
	Mat<double> p1(3, 1), p2(3, 1), p3(3, 1);
	while (Index != EOF) {
		fgets(&Index, 1, fin);
		switch(Index) {
		case 'f':
			fread(p1.data, 1, 3 * sizeof(double), fin);
			fread(p2.data, 1, 3 * sizeof(double), fin);
			fread(p3.data, 1, 3 * sizeof(double), fin);
			drawTriangle(p1, p2, p3, material);
			break;
		case 's':
			fread(p1.data, 1, sizeof(double), fin);
			fread(&R, 1, sizeof(double), fin);
			drawSphere(p1, R, material);
		default: break;
		}
	}
	fclose(fin);
}
/*--------------------------------[ 渲染 ]--------------------------------
*	[过程]:
		[1] 计算屏幕矢量、屏幕X,Y向轴
		[2] 对屏幕每个像素遍历
			[3] 计算像素矢量、光线矢量、光线追踪起点
			[4] 光线追踪算法
			[5] 基于结果绘制该像素色彩
-------------------------------------------------------------------------*/
void RayTracing::paint(const char* fileName, int sampleSt) {
	//[1]
	Mat<double> ScreenVec, ScreenXVec, ScreenYVec(3, 1);
	ScreenVec.sub(gCenter, Eye);																	//屏幕轴由眼指向屏幕中心
	ScreenYVec.getData(ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0], 1, 0).normalized();	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalized();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	//[2]
	double minDistance = 0, RayFaceDistance;
	Mat<double> PixYVec, PixXVec, PixVec, Ray, RaySt, color(3, 1);
	for (int sample = sampleSt; sample < SamplesNum; sample++) {
		writeImg(fileName); clock_t start = clock();
		Screen *= (double)sample / (sample + 1);
		for (int x = 0; x < Screen.rows; x++) {
			for (int y = 0; y < Screen.cols; y++) {
				PixVec.add(
					PixXVec.mult(x + rand() / double(RAND_MAX) - Screen.rows / 2 - 0.5, ScreenXVec),
					PixYVec.mult(y + rand() / double(RAND_MAX) - Screen.cols / 2 - 0.5, ScreenYVec)
				);//[3]
				traceRay(RaySt.add(gCenter, PixVec), Ray.add(ScreenVec, PixVec).normalized(), color.zero(), 0);	//[4][5]
				setPix(x, y, Screen(x, y) += (color *= 1.0 / (sample + 1)));
			} 
		} printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));
	}
}
/*--------------------------------[ 追踪光线 ]--------------------------------
*	[算法]:
*	[过程]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
&	[注]:distance > 1而不是> 0，是因为反射光线在接触面的精度内，来回碰自己....
-----------------------------------------------------------------------------*/
Mat<double>& RayTracing::traceRay(Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& color, int level) {
	double minDistance = DBL_MAX;
	Mat<double> intersection, intersectionTmp, FaceVec, tmp;
	Material* intersectMaterial = NULL;
	//[1][2][3]
	for (int i = 0; i < TriangleSet.size(); i++) {
		double distance = seekIntersection(TriangleSet[i], RaySt, Ray, tmp, intersectionTmp, minDistance);
		if (distance > eps) {
			minDistance = distance; intersectMaterial = TriangleSet[i].material;
			FaceVec = tmp; intersection = intersectionTmp;
		}
	}
	//[4]
	if (minDistance == DBL_MAX || intersectMaterial == NULL) return color.zero();		//Miss intersect
	if (intersectMaterial->rediateRate != 0) return color = intersectMaterial->color;	//Light Source
	if (level > maxRayLevel && rand() / double(RAND_MAX) > maxRayLevelProbability) return color.zero(); 		//Max Ray Level
	Mat<double> RayTmp;
	if (intersectMaterial->quickReflect != 0) {						//Quick Reflect: 若该点处的表面是(快速)散射面，计算点光源直接照射该点产生的颜色
		double lightCos = 0;
		FaceVec *= FaceVec.dot(Ray) > 0 ? -1 : 1;
		for (int i = 0; i < PointLight.size(); i++) {
			RayTmp.sub(PointLight[i], intersection);
			double t = (FaceVec.dot(RayTmp) / RayTmp.norm() + 1) / 2;
			lightCos = t > lightCos ? t : lightCos;
		}
		color.mult(lightCos, color.ones(3, 1));
	}
	else if (intersectMaterial->diffuseReflect != 0) {				//Diffuse Reflect
		traceRay(intersection, diffuseReflect(Ray, FaceVec, RayTmp), color.zero(), level + 1);
		color *= intersectMaterial->reflectRate;
	}
	else if (intersectMaterial->refractRate != 0) {					//Refract
		double t = refractRateBuffer; refractRateBuffer = refractRateBuffer == intersectMaterial->refractRate ? 1 : intersectMaterial->refractRate;
		refract(Ray, FaceVec, RayTmp, t, refractRateBuffer);
		traceRay(intersectionTmp.add(intersectionTmp.mult(eps, RayTmp), intersection), RayTmp, color.zero(), level + 1);
		refractRateBuffer = t;
	}
	else if (intersectMaterial->reflectRate != 0) {					//Reflect
		traceRay(intersection, reflect(Ray, FaceVec, RayTmp), color.zero(), level + 1);
		color *= intersectMaterial->reflectRate;
	}
	return color.elementMult(tmp.mult(level > maxRayLevel ? 1 / maxRayLevelProbability : 1, intersectMaterial->color));
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
double RayTracing::seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, Mat<double>& intersection, double minDistance) {
	// Sphere Seek Intersection
	if (triangle.p[2][0] == DBL_MAX) {
		// 计算ΔX、Δ
		Mat<double> RayStCenter; RayStCenter.sub(RaySt, triangle.p[0]);
		double R = triangle.p[1][0], A = Ray.dot(Ray), B = 2 * Ray.dot(RayStCenter);
		double Delta = B * B - 4 * A * (RayStCenter.dot(RayStCenter) - R * R);
		if (Delta < 0) return -1;									//有无交点
		Delta = sqrt(Delta);
		double RayFaceDistance = (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
		if (RayFaceDistance <= eps || RayFaceDistance >= minDistance) return -1;
		intersection.add(intersection.mult(RayFaceDistance, Ray), RaySt);
		FaceVec.sub(intersection, triangle.p[0]).normalized();
		return RayFaceDistance;
	}
	//[1][2]
	Mat<double> edge[2], tmp;
	FaceVec.crossProduct_(
		edge[0].sub(triangle.p[1], triangle.p[0]),
		edge[1].sub(triangle.p[2], triangle.p[0])
	).normalized();
	double t = FaceVec.dot(Ray);
	if (t == 0) return -1;											//光线与面是否平行
	double RayFaceDistance = FaceVec.dot(tmp.sub(triangle.p[0], RaySt)) / t;
	if (RayFaceDistance <= eps || RayFaceDistance >= minDistance) return -1;
	intersection.add(intersection.mult(RayFaceDistance, Ray), RaySt);
	//[3]
	tmp.sub(intersection, triangle.p[0]);
	double Dot00 = edge[0].dot(edge[0]),
		   Dot01 = edge[0].dot(edge[1]),
		   Dot11 = edge[1].dot(edge[1]),
		   Dot02 = edge[0].dot(tmp),
		   Dot12 = edge[1].dot(tmp);
	t = Dot00 * Dot11 - Dot01 * Dot01;
	double u = (Dot11 * Dot02 - Dot01 * Dot12) / t;
	double v = (Dot00 * Dot12 - Dot01 * Dot02) / t;
	return (u < 0 || u > 1 || v < 0 || v > 1 || u + v > 1) ? -1 : RayFaceDistance;
}
/*--------------------------------[ 反射 ]--------------------------------
*	[反射定律]: 入射角 == 反射角
			设面矢F, 入射光L
			Lf = L - F·2 cos<L,F>
-------------------------------------------------------------------------*/
Mat<double>& RayTracing::reflect(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& reflectRay) {
	return reflectRay.add(reflectRay.mult(-2 * faceVec.dot(incidentRay), faceVec), incidentRay).normalized();
}
/*--------------------------------[ 折射 ]--------------------------------
*	[折射定律]: n1·sinθ1 = n2·sinθ2
			设面矢F, 入射光L
			=> sin α =  n·sin <L,F>  ,  n = n入 / n折
			Lz = L + F·sin(α - <L,F>) / sinα
			   = L + F·( cos<L,F> - sqrt( 1/n² - 1 + cos²<L,F> ) )
-------------------------------------------------------------------------*/
Mat<double>& RayTracing::refract(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& refractRay, double rateIn, double rateOut) {
	double refractRate = rateIn / rateOut;
	double CosIn = faceVec.dot(incidentRay);
	double CosOut = 1 - refractRate * refractRate * (1 - CosIn * CosIn);
	if (CosOut >= 0) reflect(incidentRay, faceVec, refractRay);			//全反射
	Mat<double> tmp;
	CosOut = sqrt(CosOut);
	return refractRay.add(refractRay.mult(refractRate, incidentRay), tmp.mult((CosIn > 0 ? 1 : -1)* CosOut - refractRate * CosIn, faceVec)).normalized();
}
/*--------------------------------[ 漫反射 ]--------------------------------
*	[算法]: MonteCarlo: 随机持续采样
*	[漫反射]: 在面矢半球内，面积均匀的随机取一射线，作为反射光线.
---------------------------------------------------------------------------*/
Mat<double>& RayTracing::diffuseReflect(Mat<double>& incidentRay, Mat<double>& faceVec, Mat<double>& reflectRay) {
	double r1 = 2 * PI * rand() / double(RAND_MAX), r2 = rand() / double(RAND_MAX);
	Mat<double> tmp1(3, 1), tmp2(3, 1);
	faceVec *= faceVec.dot(incidentRay) > 0 ? -1 : 1;
	tmp1[0] = fabs(faceVec[0]) > 0.1 ? 0 : 1; tmp1[1] = tmp1[0] == 0 ? 1 : 0;
	tmp1.mult(cos(r1) * sqrt(r2), tmp1.crossProduct(tmp1, faceVec).normalized());
	tmp2.mult(sin(r1) * sqrt(r2), tmp2.crossProduct_(faceVec, tmp1));
	reflectRay.add(reflectRay.add(reflectRay.mult(sqrt(1 - r2), faceVec), tmp1), tmp2);
	return reflectRay.normalized();
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
	Mat<double> direction, rotateAxis, rotateMat(4), zAxis(3, 1), tmp; zAxis.getData(0, 0, 1);
	direction.sub(ed, st);
	if (direction[0] != 0 || direction[1] != 0) {
		GraphicsND::rotate(
			rotateAxis.crossProduct(direction, zAxis),
			-acos(tmp.dot(direction, zAxis) / direction.norm()),
			tmp.zero(3, 1), rotateMat
		);
		rotateMat.block(1, 3, 1, 3, rotateMat);
	}
	else rotateMat.E(3);
	// 画圆台
	Mat<double> stPoint, edPoint, preStPoint, preEdPoint, deltaVector(3, 1);
	for (int i = 0; i <= 360 / delta; i++) {
		deltaVector.getData(cos(i * delta * 2.0 * PI / 360), sin(i * delta * 2.0 * PI / 360), 0);
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
	p1[0] = p1[1] = p1[2] = r; p2[0] = DBL_MAX;
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