#include "RayTracing.h"

using namespace ObjectLib;

/*#############################################################################

*						光线追踪  Ray Tracing

##############################################################################*/

std::vector<Mat<>>  RayTracing::PointLight;
int RayTracing::maxRayLevel = 6;


bool RayTracing::haze = 1;
Mat<> RayTracing::haze_A{ 3 };
double RayTracing::haze_beta = 1e-3;

/*--------------------------------[ 渲染 ]--------------------------------
*	[过程]:
		[1] 计算屏幕矢量、屏幕X,Y向轴
		[2] 对屏幕每个像素遍历
			[3] 计算像素矢量、光线矢量、光线追踪起点
			[4] 光线追踪算法
			[5] 基于结果绘制该像素色彩
-------------------------------------------------------------------------*/

static void RayTracing::traceRay_func(Mat<>* ScreenXVec, Mat<>* ScreenYVec, Mat<>* center, Mat<>* direct, 
	ObjectTree* objTree, double rate, 
	Mat<>* R, Mat<>* G, Mat<>* B, int st, int ed
) {
	Mat<> SampleVec, SampleXVec, SampleYVec, Ray, RaySt, color(3);

	for (int i = st; i < ed; i++) {
		add(SampleVec,														//[3]
			mul(SampleXVec, R->i2x(i) - (*R).rows / 2.0 - 0.5, *ScreenXVec),
			mul(SampleYVec, R->i2y(i) - (*R).cols / 2.0 - 0.5, *ScreenYVec)
		);
		add(RaySt, *center, SampleVec);
		add(Ray,   *direct, SampleVec);

		RayTracing::traceRay(*objTree, RaySt, normalize(Ray), color.zero(), 0); 			//[4][5]

		mul(color, rate, color);
		(*R)(i) = (*R)(i) * (1 - rate) + color(0);
		(*G)(i) = (*G)(i) * (1 - rate) + color(1);
		(*B)(i) = (*B)(i) * (1 - rate) + color(2);
	}
}

void RayTracing::traceRay_(
	Mat<>& center, Mat<>& direct, double width, double height,
	ObjectTree& objTree,
	Mat<>& R, Mat<>& G, Mat<>& B,
	int sampleSt, int sampleEd
) {
	//[1]
	static Mat<> ScreenXVec, ScreenYVec(3);
	normalize(ScreenYVec = { -direct[1], direct[0], 0 });			//屏幕Y向轴始终与Z轴垂直,无z分量, 且与direct正交(内积为零)
	normalize(cross_(ScreenXVec, direct, ScreenYVec));				//屏幕X向轴与屏幕轴、屏幕Y向轴正交

	//[2]
	std::thread threads[4];
	int kk = R.size() / 4;

	clock_t start = clock();

	for (int sample = sampleSt; sample < sampleEd; sample++) {
		double rate = 1.0 / (sample + 1);

		for (int i = 0; i < 4; i++) {
			threads[i] = std::thread(
				traceRay_func, &ScreenXVec, &ScreenYVec, &center, &direct, &objTree, rate, &R, &G, &B, i * kk, (i + 1) * kk
			);
		}

		for (auto& thread : threads)
			thread.join();

		if (sample % 100 == 0) {
			printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));
			start = clock();
		}
	}
}


void RayTracing::traceRay(
	Mat<>& center, Mat<>& direct, double preSize,
	ObjectTree& objTree,
	Mat<>& R, Mat<>& G, Mat<>& B,
	int sampleSt, int sampleEd
) {
	//[1]
	static Mat<> ScreenXVec, ScreenYVec(3);
	normalize(ScreenYVec = { -direct[1], direct[0], 0 });			//屏幕Y向轴始终与Z轴垂直,无z分量, 且与direct正交(内积为零)
	normalize(cross_(ScreenXVec, direct, ScreenYVec));				//屏幕X向轴与屏幕轴、屏幕Y向轴正交

	//[2]
	static Mat<> SampleYVec, SampleXVec, SampleVec, Ray, RaySt, color(3);
	clock_t start = clock();

	for (int sample = sampleSt; sample < sampleEd; sample++) {
		double rate = 1.0 / (sample + 1),
			randX = RAND_DBL, 
			randY = RAND_DBL;

		for (int x = 0; x < R.rows; x++) {
			for (int y = 0; y < R.cols; y++) {
				add(SampleVec,														//[3]
					mul(SampleXVec, (x + randX - R.rows / 2.0 - 0.5) * preSize, ScreenXVec),
					mul(SampleYVec, (y + randY - R.cols / 2.0 - 0.5) * preSize, ScreenYVec)
				); 
				add(RaySt, center, SampleVec);
				add(Ray,   direct, SampleVec);

				traceRay(objTree, RaySt, normalize(Ray), color.zero(), 0); 			//[4][5]
				
				mul(color, rate, color);
				R(x, y) = R(x, y) * (1 - rate) + color(0);
				G(x, y) = G(x, y) * (1 - rate) + color(1);
				B(x, y) = B(x, y) * (1 - rate) + color(2);
			} 
		}

		if (sample % 10 == 0) {
			printf("%d\ttime:%f sec\n", sample, (clock() - start) / double(CLK_TCK));
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
Mat<>& RayTracing::traceRay(ObjectTree& objTree, Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level) 
{
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
		static Mat<> faceVec(3), tmp;

		add(RaySt, RaySt, mul(tmp, dis, Ray));
		FaceVector(*obj, RaySt, faceVec);

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
				normalize(sub(tmp, PointLight[i], RaySt));
				//t = dot(faceVec, tmp);
				t = (dot(faceVec, tmp) + 1) / 2;
				lightCos = t > lightCos ? t : lightCos;
			} 
			mul(color, lightCos, color = 1);
		}

		//漫反射
		else if (material->diffuseReflect) {
			diffuseReflect(Ray0, faceVec, Ray);
			traceRay(objTree, RaySt, Ray, color, level + 1);
			mul(color, material->reflectLossRate, color);
		}

		//反射
		else if (RAND_DBL < material->reflect) {
			reflect(Ray0, faceVec, Ray);
			traceRay(objTree, RaySt, Ray, color, level + 1);
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
			traceRay(objTree, add(RaySt, RaySt, mul(tmp, EPS, Ray)), Ray, color, level + 1);
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
	//if(haze) 
		//Haze(color, color, haze_A = 1, dis, haze_beta);

	return color;
}
