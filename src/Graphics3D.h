#ifndef GRAPHICS3D_H
#define GRAPHICS3D_H
#include "Graphics.h"
struct GraphicsCoor {
	INT32S x, y, z;
};
class Graphics3D
{
public:
	/*---------------- 基础参数 ----------------*/
	Graphics* g;															//核心图形学类
	GraphicsCoor gSize;														//图大小
	RGB PaintColor;															//画笔颜色
	INT32S PaintSize = 0;													//画笔大小
	Mat<FP64> gM;															//变换矩阵
	/*---------------- 底层 ----------------*/
	void init();															//初始化
	void clear();															//清屏
	/*---------------- DRAW ----------------*/
	void drawPoint(GraphicsCoor p0);										//画点
	void drawLine(GraphicsCoor sp, GraphicsCoor ep);						//画直线
	void drawPlane(GraphicsCoor p1, GraphicsCoor p2, GraphicsCoor p3);		//画平面
	void drawSphere(GraphicsCoor p0, INT32S r);								//画球
	void drawEllipsoid(GraphicsCoor p0, GraphicsCoor r);					//画椭球
	void drawFace();														//画曲面
	void drawBezier();														//画贝塞尔曲面
	void drawCopy(GraphicsCoor p0, RGB* gt, GraphicsCoor size);				//复制别的3D图
	void fill(GraphicsCoor sp, GraphicsCoor ep, RGB color);					//填充
	void floodfill(GraphicsCoor p0, RGB color);								//泛滥填充
	void drawChar(GraphicsCoor p0, CHAR charac);							//显示字符
	void drawString(GraphicsCoor p0, const CHAR* str, INT32S n);			//显示字符串
	void drawNum(GraphicsCoor p0, FP64 num);								//显示数字
	/*---------------- SET ----------------*/
	BOOL judgeOutRange(INT32S x0, INT32S y0);								//判断坐标是否过界
	void setSize();															//设置窗口尺寸
};

#endif