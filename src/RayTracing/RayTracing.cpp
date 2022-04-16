#include "RayTracing.h"

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
	//[0]

	objTree.build();
	//[1]
	static Mat<> ScreenVec, ScreenXVec, ScreenYVec(3);
	sub(ScreenVec, gCenter, Eye);															//屏幕轴由眼指向屏幕中心
	ScreenYVec = { ScreenVec[0] == 0 ? 0 : -ScreenVec[1] / ScreenVec[0], 1, 0 };			//屏幕Y向轴始终与Z轴垂直,无z分量
	normalize(ScreenYVec);
	cross(ScreenXVec, ScreenVec, ScreenYVec);						//屏幕X向轴与屏幕轴、屏幕Y向轴正交
	normalize(ScreenXVec);

	//[2]
	static Mat<> PixYVec, PixXVec, PixVec, Ray, RaySt, color(3); clock_t start = clock();
	for (int sample = sampleSt; sample < sampleEd; sample++) {
		double rate = 1.0 / (sample + 1), randX = RAND_DBL, randY = RAND_DBL;

		for (int x = 0; x < Screen.rows; x++) {
			for (int y = 0; y < Screen.cols; y++) {
				add(														//[3]
					PixVec,
					mul(PixXVec, x + randX - Screen.rows / 2 - 0.5, ScreenXVec),
					mul(PixYVec, y + randY - Screen.cols / 2 - 0.5, ScreenYVec)
				); 
				traceRay(														//[4][5]
					add(RaySt, gCenter,   PixVec),
					normalize(add(Ray, ScreenVec, PixVec)),
					color.zero(), 0
				); 
				
				mul(Screen(x, y), 1 - rate, Screen(x, y));
				add(Screen(x, y), Screen(x, y), mul(color, rate, color));
			} 
		}

		if (sample % 100 == 0) {
			printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));

			for (int x = 0; x < Screen.rows; x++) 
				for (int y = 0; y < Screen.cols; y++) 
					setPix(x, y, Screen(x, y));

			GraphicsIO::ppmWrite(fileName, ScreenPix); 
			start = clock();
		}
	}
}

/******************************************************************************
*						追踪光线
*	[步骤]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
&	[注]:distance > 1而不是> 0，是因为反射光线在接触面的精度内，来回碰自己....
******************************************************************************/
Mat<>& RayTracing::traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level) {
	//搜索光线-对象的最近交点距离
	Object* obj; 
	double dis = objTree.seekIntersection(RaySt, Ray, obj); 

	//无交点
	if (dis == DBL_MAX)
		color.zero();

	//光源
	else if (obj->material->rediate) 
		color = obj->material->color;

	//最大追踪层数
	else if (level > maxRayLevel) 
		return color;

	else {
		Material* material = obj->material;

		//计算面矢、交点
		static Mat<> faceVec(3), tmp, tmp2;

		add(RaySt, RaySt, mul(tmp, dis, Ray));

		switch (obj->type) {
		case PLANE:		faceVec = *(Mat<>*)obj->v[0]; break;
		case CIRCLE:	faceVec = *(Mat<>*)obj->v[1]; break;
		case TRIANGLE:	
			normalize(cross_(
				faceVec,
				sub(tmp, *(Mat<>*)obj->v[1], *(Mat<>*)obj->v[0]),
				sub(tmp2, *(Mat<>*)obj->v[2], *(Mat<>*)obj->v[0])
			)); break;
		case PLANESHAPE:faceVec = *(Mat<>*)obj->v[1]; break;
		case SPHERE:	normalize(sub(faceVec, RaySt, *(Mat<>*)obj->v[0])); break;
		case CUBOID:	 
			if      (fabs(RaySt[0] - (*(Mat<>*)obj->v[0])[0]) < EPS || fabs(RaySt[0] - (*(Mat<>*)obj->v[1])[0]) < EPS) faceVec = { 1, 0, 0 };
			else if (fabs(RaySt[1] - (*(Mat<>*)obj->v[0])[1]) < EPS || fabs(RaySt[1] - (*(Mat<>*)obj->v[1])[1]) < EPS) faceVec = { 0, 1, 0 };
			else if (fabs(RaySt[2] - (*(Mat<>*)obj->v[0])[2]) < EPS || fabs(RaySt[2] - (*(Mat<>*)obj->v[1])[2]) < EPS) faceVec = { 0, 0, 1 };
			break;
		}

		//色散补丁
		static int refractColorIndex; 
		static double refractRateBuf; 
		static bool isChromaticDisperson;
		if (level == 0) 
			refractColorIndex = RAND_DBL * 3, refractRateBuf = 1, isChromaticDisperson = 0;
		
		static Mat<> Ray0; 
		Ray0 = Ray;

		//快速反射
		if (material->quickReflect) {								//Reflect Quick: 计算点光源直接照射该点产生的颜色
			double lightCos = 0, t;
			mul(faceVec, dot(faceVec, Ray) > 0 ? -1 : 1, faceVec);
			for (int i = 0; i < PointLight.size(); i++) {
				t = dot(faceVec, normalize(sub(tmp, PointLight[i], RaySt)));
				lightCos = t > lightCos ? t : lightCos;
			} 
			mul(color, lightCos, color = 1);
		}

		//漫反射
		else if (material->diffuseReflect) {
			diffuseReflect(Ray0, faceVec, Ray);
			traceRay(RaySt, Ray, color, level + 1);
			mul(color, material->reflectLossRate, color);
		}

		//反射
		else if (RAND_DBL < material->reflect) {
			reflect(Ray0, faceVec, Ray);
			traceRay(RaySt, Ray, color, level + 1);
			mul(color, material->reflectLossRate, color);
		}

		//折射
		else {
			//折射率补丁
			if (material->refractRate[0] != material->refractRate[1] 
			 || material->refractRate[0] != material->refractRate[2]) 
				isChromaticDisperson = 1;
			double t = refractRateBuf; 
			refractRateBuf = refractRateBuf == material->refractRate[refractColorIndex] ? 
				1 : material->refractRate[refractColorIndex];

			refract(Ray0, faceVec, Ray, t, refractRateBuf);
			traceRay(add(RaySt, RaySt, mul(tmp, EPS, Ray)), Ray, color, level + 1);
			mul(color, material->refractLossRate, color);
		}
		//色散补丁
		if (level == 0 && isChromaticDisperson) { 
			double t = color[refractColorIndex]; 
			color.zero();
			color[refractColorIndex] = 3 * t;
		}

		elementMul(color, color, obj->material->color);
	}

	//雾
	if(haze) 
		Haze(color, color, haze_A, dis, haze_beta);

	return color;
}
