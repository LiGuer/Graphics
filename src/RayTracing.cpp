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
/*#############################################################################

*						几何光学  Geometrical Optics

-------------------------------------------------------------------------------
*	[基本定律]:
		[1] 光沿直线传播.
		[2] 两束光传播中互不干扰,会聚于同一点时光强简单相加.
		[3] 介质分界面光传播:
			[3.1] 反射. 入射角 == 反射角
			[3.2] 折射. n1·sinθ1 = n2·sinθ2
			[3.3] 漫反射
		[4] 光路可逆
##############################################################################*/
/******************************************************************************
*						反射
*	[反射定律]: 入射角 == 反射角
			设面矢F, 入射光L
			Lf = L - F·2 cos<L,F>
******************************************************************************/
Mat<>& reflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
	return RayO.add(RayO.mul(-2 * faceVec.dot(RayI), faceVec), RayI).normalize();
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
Mat<>& refract(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO) {
	double k = rateI / rateO,
		CosI = faceVec.dot(RayI),
		CosO = 1 - pow(k, 2) * (1 - pow(CosI, 2));
	return CosO < 0 ? reflect(RayI, faceVec, RayO) :				//全反射
		RayO.add(RayI, RayO.mul(-CosI - (CosI > 0 ? -1 : 1)* sqrt(CosO) / k, faceVec)).normalize();
}
/******************************************************************************
*						漫反射
*	[算法]: MonteCarlo: 随机持续采样
*	[漫反射]: 在面矢半球内，面积均匀的随机取一射线，作为反射光线.
******************************************************************************/
Mat<>& diffuseReflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
	double r1 = 2 * PI * RAND_DBL, r2 = RAND_DBL;
	static Mat<> t(3), u, v;
	faceVec *= faceVec.dot(RayI) > 0 ? -1 : 1;
	t[0] = fabs(faceVec[0]) > 0.1 ? 0 : 1;
	t[1] = t[0] == 0 ? 1 : 0;
	u.mul(cos(r1) * sqrt(r2), u.cross_(t, faceVec).normalize());
	v.mul(sin(r1) * sqrt(r2), v.cross_(faceVec, u).normalize());
	return RayO.add(RayO.mul(sqrt(1 - r2), faceVec), u += v).normalize();
}

/******************************************************************************
						求交点
*	[算法]:
-------------------------------------------------------------------------------
		[射线平面交点]
		平面方程: Ax + By + Cz + D = 0
		直线方程: (x - x0) / a = (y - y0) / b = (z - z0) / c = K
			联立平面,直线方程, 解K, K即直线到射线平面交点的距离
			x = K a + x0 , ...
			A(K a + x0) + B(K b + y0) + C(K c + z0) + D = 0
			=> K =  - (A x0 + B y0 + C z0) / (A a + B b + C c)
-------------------------------------------------------------------------------
		[射线球面交点]
		球方程: (X - Xs)² + (Y - Ys)² + (Z - Zs)² = R²
		线球交点: K²(Al² + Bl² + Cl²) + 2K(Al ΔX + Bl ΔY + Cl ΔZ) + (ΔX² + ΔY² + ΔZ² - R²) = 0
				  ΔX = Xl - Xs
				  Δ = b² - 4ac = 4(Al ΔX + Bl ΔY + Cl ΔZ)² - 4(Al² + Bl² + Cl²)(ΔX² + ΔY² + ΔZ² - R²)
				  若Δ≥0 有交点.
				  K = ( -b ± sqrt(Δ) ) / 2a	即光线走过线距离
-------------------------------------------------------------------------------
		[射线三角形交点] Moller-Trumbore方法(1997)
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
-------------------------------------------------------------------------------
		[射线长方体交点]
******************************************************************************/
double RayPlane(Mat<>& RaySt, Mat<>& Ray, double& A, double& B, double& C, double& D) {
	double t = A * Ray[0] + B * Ray[1] + C * Ray[2];
	if (t == 0) return DBL_MAX;
	double d = -(A * RaySt[0] + B * RaySt[1] + C * RaySt[2] + D) / t;
	return d > 0 ? d : DBL_MAX;
}
double RayCircle(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, double& R, Mat<>& normal) {
	double
		D = -(normal[0] * Center[0] + normal[1] * Center[1] + normal[2] * Center[2]),
		d = RayPlane(RaySt, Ray, normal[0], normal[1], normal[2], D);
	static Mat<> tmp; tmp.add(RaySt, tmp.mul(d, Ray)); 
	return (tmp -= Center).norm() <= R ? d : DBL_MAX;
}
double RayTriangle(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	static Mat<> edge[2], tmp, p, q;
	edge[0].sub(p2, p1);
	edge[1].sub(p3, p1);
	// p & a & tmp
	static double a, u, v;
	a = p.cross_(Ray, edge[1]).dot(edge[0]);
	if (a > 0) tmp.sub(RaySt, p1);
	else       tmp.sub(p1, RaySt), a = -a;
	if (a < 1e-4)		return DBL_MAX;								//射线与三角面平行
	// u & q & v
	u = p.dot(tmp) / a;
	if (u < 0 || u > 1)	return DBL_MAX;
	v = q.cross_(tmp, edge[0]).dot(Ray) / a;
	return (v < 0 || u + v > 1) ? DBL_MAX : q.dot(edge[1]) / a;
}
double RaySphere(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R) {
	static Mat<> RayStCenter; RayStCenter.sub(RaySt, center);
	double 
		A = Ray.dot(Ray),
		B = 2 * Ray.dot(RayStCenter),
		Delta = B * B - 4 * A * (RayStCenter.dot(RayStCenter) - R * R);
	if (Delta < 0) return DBL_MAX;									//有无交点
	Delta = sqrt(Delta);
	return (-B + (-B - Delta > 0 ? -Delta : Delta)) / (2 * A);
}
double RayCuboid(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	return DBL_MAX;
}
double RayCuboid(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax) {
	double t0 = -DBL_MAX, t1 = DBL_MAX;
	for (int dim = 0; dim < 3; dim++) {
		if (Ray[dim] == 0) {
			return DBL_MAX;
		}
		double
			t0t = (pmin[dim] - RaySt[dim]) / Ray[dim],
			t1t = (pmax[dim] - RaySt[dim]) / Ray[dim];
		if (t0t > t1t) std::swap(t0t, t1t);
		t0 = std::max(t0, t0t);
		t1 = std::min(t1, t1t);
		if (t0 > t1 || t1 < 0) return DBL_MAX;
	}
	return t0 >= 0 ? t0 : t1;
}

/*#############################################################################

*						光线追踪  Ray Tracing

##############################################################################*/
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
	ScreenVec. sub(gCenter, Eye);															//屏幕轴由眼指向屏幕中心
	ScreenYVec.set(ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0], 1, 0).normalize();	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.cross(ScreenVec, ScreenYVec).normalize();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	//[2]
	static Mat<> PixYVec, PixXVec, PixVec, Ray, RaySt, color(3); clock_t start;
	for (int sample = sampleSt; sample < sampleEd; sample++) {
		if (sample % 1 == 0) { GraphicsFileCode::ppmWrite(fileName, ScreenPix); start = clock(); }
		double rate = 1.0 / (sample + 1);
		for (int x = 0; x < Screen.rows; x++) {
			for (int y = 0; y < Screen.cols; y++) {
				PixVec.add(														//[3]
					PixXVec.mul(x + RAND_DBL - Screen.rows / 2 - 0.5, ScreenXVec),
					PixYVec.mul(y + RAND_DBL - Screen.cols / 2 - 0.5, ScreenYVec)
				);
				traceRay(														//[4][5]
					RaySt.add(gCenter,   PixVec), 
					Ray.  add(ScreenVec, PixVec).normalize(), 
					color.zero(), 0
				);
				setPix(x, y, (Screen(x, y) *= 1 - rate) += (color *= rate));
			} 
		} if (sample % 1 == 0) printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));
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
	Object* minDisOb = NULL;
	for (int i = 0; i < ObjectSet.size(); i++) {
		dis    = seekIntersection(ObjectSet[i], RaySt, Ray);
		minDis = (dis > eps && dis < minDis) ? minDisOb = &ObjectSet[i], dis : minDis;
	}
	if (minDis == DBL_MAX)			return color;					//Miss intersect
	Material* material = minDisOb->material;
	if (material->rediate != 0)		return color = material->color;	//Light Source
	if (level > maxRayLevel)		return color;					//Max Ray Level
	//[4] RaySt & FaceVec
	static Mat<> faceVec, RayTmp, tmp;
	{
		RaySt += (tmp.mul(minDis, Ray));
		switch (minDisOb->type) {
		case PLANE:		faceVec = *(Mat<>*)minDisOb->v[0]; break;
		case CIRCLE:	faceVec = *(Mat<>*)minDisOb->v[1]; break;
		case TRIANGLE:	faceVec.cross_(
							   tmp.sub(*(Mat<>*)minDisOb->v[1], *(Mat<>*)minDisOb->v[0]),
							RayTmp.sub(*(Mat<>*)minDisOb->v[2], *(Mat<>*)minDisOb->v[0])
						).normalize(); break;
		case SPHERE:	faceVec.sub(RaySt, *(Mat<>*)minDisOb->v[0]).normalize(); break;
		case CUBOID:	if (fabs(RaySt[0] - (*(Mat<>*)minDisOb->v[0])[0]) < eps || fabs(RaySt[0] - (*(Mat<>*)minDisOb->v[0])[0]) < eps) faceVec = { 1, 0, 0 };
				   else if (fabs(RaySt[1] - (*(Mat<>*)minDisOb->v[0])[1]) < eps || fabs(RaySt[1] - (*(Mat<>*)minDisOb->v[0])[1]) < eps) faceVec = { 0, 1, 0 };
				   else if (fabs(RaySt[2] - (*(Mat<>*)minDisOb->v[0])[2]) < eps || fabs(RaySt[2] - (*(Mat<>*)minDisOb->v[0])[2]) < eps) faceVec = { 0, 0, 1 };
						break;
		}
	}
	//[5]
	static int refractColorIndex; static double refractRateBuf; static bool isChromaticDisperson;
	if (level == 0) refractColorIndex = RAND_DBL * 3, refractRateBuf = 1, isChromaticDisperson = 0;
	RayTmp = Ray;
	if (material->quickReflect) {								//Reflect Quick: 计算点光源直接照射该点产生的颜色
		double lightCos = 0, t;
		faceVec *= faceVec.dot(Ray) > 0 ? -1 : 1;
		for (int i = 0; i < PointLight.size(); i++) {
			t = faceVec.dot(tmp.sub(PointLight[i], RaySt).normalize());
			lightCos = t > lightCos ? t : lightCos;
		} color.fill(1) *= lightCos;
	}
	else if (material->diffuseReflect) {						//Reflect Diffuse
		traceRay(RaySt, diffuseReflect(RayTmp, faceVec, Ray), color, level + 1);
		color *= material->reflectLossRate;
	}
	else if (RAND_DBL < material->reflect) {					//Reflect
		traceRay(RaySt, reflect       (RayTmp, faceVec, Ray), color, level + 1);
		color *= material->reflectLossRate;
	}
	else{														//Refract
		if(material->refractRate [0] != material->refractRate[1] || material->refractRate[0] != material->refractRate[2]) isChromaticDisperson = 1;
		double t = refractRateBuf; refractRateBuf = refractRateBuf == material->refractRate[refractColorIndex] ? 1 : material->refractRate[refractColorIndex];
		refract(RayTmp, faceVec, Ray, t, refractRateBuf);
		traceRay(RaySt += (tmp.mul(eps, Ray)), Ray, color, level + 1);
		color *= material->refractLossRate; 
	}
	if (level == 0 && isChromaticDisperson) { double t = color[refractColorIndex]; color.zero()[refractColorIndex] = 3 * t; }
	return color.elementMul(material->color);
}
double RayTracing::seekIntersection(Object& ob, Mat<>& RaySt, Mat<>& Ray) {
	switch (ob.type) {
	case PLANE:		return RayPlane		(RaySt, Ray,(*(Mat<>*)ob.v[0])[0], (*(Mat<>*)ob.v[0])[1], (*(Mat<>*)ob.v[0])[2], *(double*)ob.v[1]);
	case CIRCLE:	return RayCircle	(RaySt, Ray, *(Mat<>*)ob.v[0], *(double*)ob.v[2], *(Mat<>*)ob.v[1]);
	case TRIANGLE:	return RayTriangle	(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>* )ob.v[1], *(Mat<>*)ob.v[2]);
	case SPHERE:	return RaySphere	(RaySt, Ray, *(Mat<>*)ob.v[0], *(double* )ob.v[0]);
	case CUBOID:	return RayCuboid	(RaySt, Ray, *(Mat<>*)ob.v[0], *(Mat<>* )ob.v[1]);
	}
}
/*--------------------------------[ add Object ]--------------------------------*/
void RayTracing::addPlane(Mat<>& k, Mat<>& p0, Material* material) {
	Object ob; ob.type = PLANE;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = k;(*(Mat<>*)ob.v[0]).normalize();
	ob.v[1] = new double;	*(double*)ob.v[1] = -((*(Mat<>*)ob.v[0])[0] * p0[0] + (*(Mat<>*)ob.v[0])[1] * p0[1] + (*(Mat<>*)ob.v[0])[2] * p0[2]);
	ob.material = material;
	ObjectSet.push_back(ob);
}
void RayTracing::addCircle(Mat<>& center, double R, Mat<>& n, Material* material) {
	Object ob; ob.type = CIRCLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new Mat<>;	*(Mat<>*) ob.v[1] = n; (*(Mat<>*)ob.v[1]).normalize();
	ob.v[2] = new double;	*(double*)ob.v[2] = R;
	ob.material = material;
	ObjectSet.push_back(ob);
}
void RayTracing::addTriangle(Mat<>& p1, Mat<>& p2, Mat<>& p3, Material* material) {
	Object ob; ob.type = TRIANGLE;
	ob.v = (void**)calloc(3, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = p1;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = p2;
	ob.v[2] = new Mat<>;	*(Mat<>*)ob.v[2] = p3;
	ob.material = material;
	ObjectSet.push_back(ob);
}
void RayTracing::addSphere(Mat<>& center, double r, Material* material) {
	Object ob; ob.type = SPHERE;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*) ob.v[0] = center;
	ob.v[1] = new double;	*(double*)ob.v[1] = r;
	ob.material = material;
	ObjectSet.push_back(ob);
}
void RayTracing::addCuboid(Mat<>& pmin, Mat<>& pmax, Material* material){
	Object ob; ob.type = CUBOID;
	ob.v = (void**)calloc(2, sizeof(void*));
	ob.v[0] = new Mat<>;	*(Mat<>*)ob.v[0] = pmin;
	ob.v[1] = new Mat<>;	*(Mat<>*)ob.v[1] = pmax;
	ob.material = material;
	ObjectSet.push_back(ob);
}
void RayTracing::addStl(const char* file, Mat<>& center, double size, Material** material) {
	Mat<> p0(3), p1(3), p2(3), p3(3), p4(3), p5(3), p6(3); Mat<short> a;
	GraphicsFileCode::stlRead(file, p0, p1, p2, p3, a);
	for (int i = 0; i < p0.cols; i++) {
		addTriangle(
			((p4 = { p1(0,i), p1(1,i), p1(2,i) }) *= size) += center,
			((p5 = { p2(0,i), p2(1,i), p2(2,i) }) *= size) += center,
			((p6 = { p3(0,i), p3(1,i), p3(2,i) }) *= size) += center,
			material[a[i]]
		);
	}
}