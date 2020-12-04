#ifndef GRAPHICS3D_H
#define GRAPHICS3D_H
#include "Graphics.h"
class Graphics3D
{
public:
	/*---------------- 基础参数 ----------------*/
	Graphics* g = NULL;														//核心图形学类
	Mat<FP64> WindowSize{ 2,1 };											//窗口尺寸
	Mat<FP64> TransformMat;													//变换矩阵
	/*---------------- 底层 ----------------*/
	Graphics3D(int WindowSize_Width, int WindowSize_height);				//构造函数
	~Graphics3D();															//析构函数
	void init(int WindowSize_Width, int WindowSize_height);					//初始化
	void clear();															//清屏
	/*---------------- DRAW ----------------*/
	void drawPoint(Mat<double>& p0);										//画点
	void drawLine(Mat<double>& sp, Mat<double>& ep);						//画直线
	void drawPlane(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3);		//画平面
	void drawSphere(Mat<double>& p0, INT32S r);								//画球
	void drawEllipsoid(Mat<double>& p0, Mat<double>& r);					//画椭球
	void drawFace();														//画曲面
	void drawBezier();														//画贝塞尔曲面
	void fill(Mat<double>& sp, Mat<double>& ep, RGB color);					//填充
	void fillTriangle(Mat<double> p0[]);									//三角形填充
	void fillflood(Mat<double>& p0, RGB color);								//泛滥填充
	void drawChar(Mat<double>& p0, CHAR charac);							//显示字符
	void drawString(Mat<double>& p0, const CHAR* str, INT32S n);			//显示字符串
	void drawNum(Mat<double>& p0, FP64 num);								//显示数字
	/*---------------- Transformation ----------------*/
	void translation(Mat<double>& delta);									//平移
	void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center);//三维旋转-四元数
	void scaling(Mat<double>& scale, Mat<double>& center);					//缩放
	/*---------------- SET ----------------*/
	BOOL judgeOutRange(INT32S x0, INT32S y0);								//判断坐标是否过界
	void setSize();															//设置窗口尺寸
};

#endif