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
#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <queue>
#include "font.h"
#include "../LiGu_AlgorithmLib/Mat.h"
const double PI = 3.141592653598;
/******************************************************************************
*                    定义与编译器无关的数据类型
******************************************************************************/
typedef unsigned char  INT8U;			/* Unsigned  8 bit quantity       */
typedef signed   char  INT8S;			/* Signed    8 bit quantity       */
typedef unsigned short INT16U;			/* Unsigned 16 bit quantity       */
typedef signed   short INT16S;			/* Signed   16 bit quantity       */
typedef unsigned int   INT32U;			/* Unsigned 32 bit quantity       */
typedef signed int     INT32S;			/* Signed   32 bit quantity       */
typedef long long      INT64S;			/* Signed   64 bit quantity       */
typedef float          FP32;			/* Single precision floating point*/
typedef double         FP64;			/* Double precision floating point*/
typedef char    CHAR;			/* 字符       */
typedef INT32U  RGB;
typedef INT8U  RGBBASIC;
/******************************************************************************
*                    Graphics 计算机图形学
******************************************************************************/
class Graphics {
public:
	/*---------------- 基础参数 ----------------*/
	INT32S gWidth = 100, gHeight = 100;										//窗口尺寸
	RGBBASIC* Map = NULL;														//图
	RGB PaintColor = 0xFFFFFF;												//画笔颜色
	INT32S PaintSize = 0, FontSize=16;										//画笔大小//字符大小
	Mat<FP64> gM;															//变换矩阵
	/*---------------- 常数 ----------------*/
	const RGB TRANSPARENT = 0xFFFFFFFF;										//RGB:透明
	/*---------------- 底层 ----------------*/
	Graphics() { ; }
	Graphics(INT32S width, INT32S height) { init(width, height); }
	~Graphics() { free(Map); }												//析构函数
	void init();															//初始化
	void init(INT32S width, INT32S height);									//初始化
	void clear(RGB color);	 												//清屏
	void setPoint(INT32S x, INT32S y, RGB color);							//底层画点
	RGB  readPoint(INT32S x, INT32S y); 									//读点 
	void PicWrite(const CHAR* filename);									//存图
	void confirmTrans();													//确认变换
	/*---------------- DRAW ----------------*/
	void drawPoint(INT32S x0, INT32S y0);									//画点
	void drawLine(INT32S x1, INT32S y1, INT32S x2, INT32S y2);				//画线
	void drawCircle(INT32S x0, INT32S y0, INT32S r);					    //画圆
	void drawEllipse(INT32S x0, INT32S y0, INT32S rx,INT32S ry);			//画椭圆
	void drawTriangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2, INT32S x3, INT32S y3);//画三角形
	void drawRectangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		   	//画矩形
	void drawWave(INT32S x[], INT32S y[], INT32S n);						//画曲线
	void drawBezier(INT32S x[], INT32S y[], INT32S n);						//画贝塞尔曲线
	void drawGrid(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy);//画网格
	void drawCopy(INT32S x0, INT32S y0, RGBBASIC* gt, INT32S width, INT32S height);//复制别的图
	void fill(INT32S sx, INT32S sy, INT32S ex, INT32S ey, RGB color);		//填充单色
	void fillflood(INT32S x0, INT32S y0, RGB color);						//泛滥填充
	void fillPolygon(INT32S x[], INT32S y[], INT32S n);						//多边形填充
	void drawChar(INT32S x0, INT32S y0, CHAR charac);						//显示字符
	void drawString(INT32S x0, INT32S y0, const CHAR* str, INT32S n);		//显示字符串
	void drawNum(INT32S x0, INT32S y0, FP64 num);							//显示数字
	/*---------------- 二维变换 TRANSFORMATION ----------------*/
	void translation(INT32S dx, INT32S dy);									//平移
	void rotate(FP64 theta, INT32S x0, INT32S y0);							//旋转
	void scaling(FP64 sx, FP64 sy);											//缩放(>1直接完成变换)
	/*---------------- SET ----------------*/
	bool judgeOutRange(INT32S x0, INT32S y0);								//判断坐标是否过界
	void setSize(INT32S width, INT32S height);								//设置窗口尺寸
};
#endif
