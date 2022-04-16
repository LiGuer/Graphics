#ifndef RAY_TRACING_H
#define RAY_TRACING_H
#include <time.h>
#include <vector>
#include <algorithm>
#include "../RGB.h"
#include "../GraphicsIO.h"
#include "GeometricalOptics.h"
#include "Object.h"

using namespace GeometricalOptics;

#define PI 3.141592653589
#define EPS 10e-4
#define RAND_DBL (rand() / double(RAND_MAX))

/*---------------- 光线追踪 ----------------*/
class RayTracing {
public:
	//基础参数 
	Mat<> gCenter{ 3 }, Eye{ 3 };
	Mat<RGB>	ScreenPix;
	Mat<Mat<>>	Screen;
	int maxRayLevel = 6;
	double ScreenXSize, ScreenYSize;
	ObjectTree objTree;
	std::vector<Mat<>>  PointLight;													//点光源集(QuickReflect专用)

	bool haze = false; 
	Mat<> haze_A{ 3 };
	double haze_beta = 1;

	//函数
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }						//构造函数
	void init (int width, int height);												//初始化
	void setPix(int x, int y, Mat<>& color);										//画像素
	void paint(const char* fileName, int sampleSt = 0, int sampleEd = 0x7FFFFFFF);	//渲染
	Mat<>& traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);				//追踪光线
};

#endif