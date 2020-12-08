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
#include "Graphics.h"
/******************************************************************************
*                    Basic Function
******************************************************************************/
/* ---------------- INIT ---------------- */
void Graphics::init() {
	if (Map != NULL)free(Map);
	Map = (RGBBASIC*)malloc(sizeof(RGBBASIC) * 3 * gWidth * gHeight);
	clear(0);
	gM.E(3);
}
void Graphics::init(INT32S width, INT32S height) {
	gWidth = width; gHeight = height;
	init();
}
/* ---------------- CLEAR ---------------- */
void Graphics::clear(RGB color)
{
	if (color == TRANSPARENT || color == 0) {	//memset按字节处理，故只能处理高低字节相同的值
		memset(Map, color, sizeof(RGBBASIC) * 3 * gWidth * gHeight);
		return;
	}
	for (INT32S y = 0; y < gHeight; y++) {
		for (INT32S x = 0; x < gWidth; x++) {
			Map[(y * gWidth + x) * 3 + 0] = color >> 16;
			Map[(y * gWidth + x) * 3 + 1] = color >> 8;
			Map[(y * gWidth + x) * 3 + 2] = color;
		}
	}
}
/* ---------------- SET/READ POINT ---------------- 
*	AlphaBlend 算法,	8位ARGB色彩
** ---------------------------------------- */
void Graphics::setPoint(INT32S x, INT32S y,RGB color) {
	if (judgeOutRange(x, y))return;
	double alpha = (color >> 24) / 255.0;
	unsigned char R = color >> 16, G = color >> 8, B = color;
	Map[(y * gWidth + x) * 3 + 0] = alpha * Map[(y * gWidth + x) * 3 + 0] + (1 - alpha) * R;
	Map[(y * gWidth + x) * 3 + 1] = alpha * Map[(y * gWidth + x) * 3 + 1] + (1 - alpha) * G;
	Map[(y * gWidth + x) * 3 + 2] = alpha * Map[(y * gWidth + x) * 3 + 2] + (1 - alpha) * B;
}
RGB Graphics::readPoint(INT32S x, INT32S y) {
	if (judgeOutRange(x, y))return TRANSPARENT;
	RGB R = Map[(y * gWidth + x) * 3 + 0];
	RGB G = Map[(y * gWidth + x) * 3 + 1];
	RGB B = Map[(y * gWidth + x) * 3 + 2];
	return R * 0x10000 + G * 0x100 + B;
}
/* ---------------- PicWrite ---------------- */
void Graphics::PicWrite(const CHAR* filename) {		// 太慢
	FILE* fp = fopen(filename, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", gWidth, gHeight);// 写图片格式、宽高、最大像素值
	fwrite(Map, 1, gHeight * gWidth * 3, fp);// 写RGB数据
	fclose(fp);
}
/*---------------- CONFIRM Trans ----------------*/
void Graphics::confirmTrans()
{
	Graphics Maptemp;
	Maptemp.setSize(gWidth, gHeight);
	Maptemp.init();
	Maptemp.clear(TRANSPARENT);
	Mat<FP64> p;
	p.zero(3, 1);
	p[2] = 1;
	for (INT32S y = 0; y < gHeight; y++) {
		for (INT32S x = 0; x < gWidth; x++) {
			p[0] = (FP64)x; p[1] = (FP64)y;
			p.mult(gM, p);
			RGB t = readPoint(x, y);
			if (t != TRANSPARENT)Maptemp.setPoint(p[0], p[1], t);
			if (t != TRANSPARENT)Maptemp.setPoint(p[0] + 0.4, p[1] + 0.4, t);	//#补丁
		}
	}
	gM.E(3);
	free(Map);
	Map = Maptemp.Map;	Maptemp.Map = NULL;
}
/******************************************************************************
*                    底层无关
******************************************************************************/
/* ---------------- DRAW POINT ---------------- */
void Graphics::drawPoint(INT32S x0, INT32S y0) {
	if (judgeOutRange(x0, y0))return;
	setPoint(x0, y0, PaintColor);						//基础点(点粗==0)
	/*------ 点粗>0时 ------*/
	if (PaintSize > 0) {
		INT32S x = 0, y = PaintSize, p = 3 - (PaintSize << 1);		//初始点:天顶(0,r)//p:决策参数(r右移即乘2)
		INT32S x_step[] = { 1,1,-1,-1 }, y_step[] = { 1,-1,1,-1 };		//上下左右对称四个点
		/*------ 绘制圆 (x=0始,y=x终) ------*/
		while (x <= y) {
			for (int i = x0 - x; i <= x0 + x; i++) { setPoint(i, y0 - y, PaintColor); setPoint(i, y0 + y, PaintColor); }	//填充圆内
			for (int i = x0 - y; i <= x0 + y; i++) { setPoint(i, y0 - x, PaintColor); setPoint(i, y0 + x, PaintColor); }
			x++;
			INT32S dp = 4 * x + 6;
			if (p < 0)p += dp;
			else {														//p过临界值，y加一
				p += dp - 4 * y + 4;
				y--;
			}
		}
	}
}
/* ---------------- DRAW LINE ---------------- 
*	Bresenham Algorithm		
*	优化：
		1. 化[浮点运算]为[整数运算]：err
		2. 各方向均可绘制
** ---------------------------------------- */
void Graphics::drawLine(INT32S x1, INT32S y1, INT32S x2, INT32S y2) {
	INT32S err[2] = { 0 }, inc[2] = { 0 };
	INT32S delta[2] = { x2 - x1, y2 - y1 };						//计算坐标增量
	INT32S x = x1, y = y1;
	//设置x单步方向	
	if (delta[0] > 0)inc[0] = 1; 								//向右
	else if (delta[0] == 0)inc[0] = 0;							//垂直
	else { inc[0] = -1; delta[0] = -delta[0]; }					//向左
	//设置y单步方向	
	if (delta[1] > 0)inc[1] = 1;								//向上
	else if (delta[1] == 0)inc[1] = 0;							//水平
	else { inc[1] = -1; delta[1] = -delta[1]; }					//向下

	INT32S distance = delta[0] > delta[1] ? delta[0] : delta[1];//总步数
	//画线
	for (INT32S i = 0; i <= distance + 1; i++) {				
		drawPoint(x, y);										//唯一输出：画点
		err[0] += delta[0];
		err[1] += delta[1];
		if (err[0] > distance) {
			err[0] -= distance;
			x += inc[0];										//x走一步
		}
		if (err[1] > distance) {
			err[1] -= distance;
			y += inc[1];										//y走一步
		}
	}
}
/* ---------------- DRAW TRIANGLE ---------------- */
void Graphics::drawPolygon(INT32S x[], INT32S y[], INT32S n)
{
	for (int i = 0; i < n; i++) 
		drawLine(x[i], y[i], x[(i + 1) % n], y[(i + 1) % n]);
}
/* ---------------- DRAW RACTANGLE ---------------- */
//画矩形	  
//(x1,y1),(x2,y2):矩形的对角坐标
void Graphics::drawRectangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2)
{
	drawLine(x1, y1, x2, y1);
	drawLine(x1, y1, x1, y2);
	drawLine(x1, y2, x2, y2);
	drawLine(x2, y1, x2, y2);
}
/* ---------------- DRAW CIRCLE ---------------
*	Bresenham Algorithm 画圆
*	充分利用圆的对称性: 1/8圆 (主: 45`-90`,斜率<=1)
*
*	通过设定在每一步取样步骤中寻找最接近圆周像素的决策参数，
*	可以将Bresenham线法移植为画圆算法。
** ---------------------------------------- */
void Graphics::drawCircle(INT32S x0, INT32S y0, INT32S r)
{
	INT32S x = 0, y = r, p = 3 - (r << 1);		//初始点:天顶(0,r)//p:决策参数(r右移即乘2)
	INT32S x_step[] = { 1,1,-1,-1 }, y_step[] = { 1,-1,1,-1 };		//上下左右对称四个点
	/*------ 绘制圆 (x=0始,y=x终) ------*/
	while (x <= y){
		for (int i = 0; i < 4; i++) {								//八分圆对称绘制
			drawPoint(x0 + x * x_step[i], y0 + y * y_step[i]);
			drawPoint(x0 + y * x_step[i], y0 + x * y_step[i]);
		}
		x++;
		INT32S dp = 4 * x + 6;
		if (p < 0)p += dp;
		else{														//p过临界值，y加一
			p += dp - 4 * y + 4;
			y--;
		}
	}
}
/* ---------------- DRAW ELLIPSE ---------------
*	中点算法
*	椭圆的对称性: 1/4椭圆
*	区域1:斜率<=1	区域2:斜率>1
*	椭圆函数:f = x²ry² + y²rx² - rx²ry² = 0	(>0界外，<0界内)
*	椭圆斜率:   dy       2 ry² * x
*	           ---- = - ------------
*	            dx       2 rx² * y
*	
*	决策参数p:侯选位置(x+1，y)与(x+1，y-1)中点，带入f所得的值。
*	由[决策参数]确定，下一步是(x+1，y)or(x+1，y-1),
*	若p < 0, 则中点在界内，(x+1，y)更合适; 否则另一个。
** ---------------------------------------- */
void Graphics::drawEllipse(INT32S x0, INT32S y0, INT32S rx, INT32S ry)
{
	INT64S rx2 = rx * rx, ry2 = ry * ry;
	INT32S x = 0, y = ry;											//初始点:天顶
	INT64S p = ry2 + rx2 * (0.25 - ry);								//p:决策参数
	INT32S x_step[] = { 1,1,-1,-1 }, y_step[] = { 1,-1,1,-1 };		//上下左右对称四个点
	/*------ 绘制椭圆 (区域1) ------*/
	while (ry2 * x < rx2 * y) {
		for (int i = 0; i < 4; i++) {								//四分对称绘制
			drawPoint(x0 + x * x_step[i], y0 + y * y_step[i]);
		}
		x++;
		INT64S dp = (1 + 2 * x) * ry2;
		if (p < 0)p += dp;
		else {														//p过临界值，y加一
			y--;
			p += dp - 2 * rx2 * y;
		}
	}
	/*------ 绘制椭圆 (区域2) ------*/
	p = ry2 * (x + 0.5) * (x + 0.5) + rx2 * (y - 1) * (y - 1) - ry2 * rx2;
	while (y >= 0) {
		for (int i = 0; i < 4; i++) {								//八分圆对称绘制
			drawPoint(x0 + x * x_step[i], y0 + y * y_step[i]);
		}
		y--;
		INT64S dp = (1 - 2 * y) * rx2;
		if (p > 0)p += dp;
		else {														//p过临界值，y加一
			x++;
			p += dp + 2 * ry2 * x;
		}
	}
}
/* ---------------- DRAW GRID ---------------- */
void Graphics::drawGrid(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy)
{
	if (dx > 0) {
		for (int x = sx; x <= ex; x += dx) {
			if (judgeOutRange(x, 0))continue;
			drawLine(x, sy, x, ey);
		}
	}
	else {
		for (int x = sx; x >= ex; x += dx) {
			if (judgeOutRange(x, 0))continue;
			drawLine(x, sy, x, ey);
		}
	}
	if (dy > 0) {
		for (int y = sy; y <= ey; y += dy) {
			if (judgeOutRange(0, y))continue;
			drawLine(sx, y, ex, y);
		}
	}
	else {
		for (int y = sy; y >= ey; y += dy) {
			if (judgeOutRange(0, y))continue;
			drawLine(sx, y, ex, y);
		}
	}
}
/* ---------------- DRAW WAVE ---------------- */
void Graphics::drawWave(INT32S x[], INT32S y[], INT32S n) {
	INT32S xt, yt;
	for (INT32S i = 0; i < n; i++) {
		if (i != 0)drawLine(xt, yt, x[i], y[i]);
		xt = x[i], yt = y[i];	
	}
}
/* ---------------- DRAW BEZIER CURVE ---------------- */
void Graphics::drawBezier(INT32S xCtrl[], INT32S yCtrl[], INT32S n)
{
	INT32S N = gWidth+ gHeight;					//#待优化
	FP64 C[50];
	for (INT32S i = 0; i < n; i++) {
		setPoint(xCtrl[i], yCtrl[i],0xFFFFFF);
	}
	for (INT32S i = 0; i < n; i++) {
		C[i] = 1;
		for (INT32S j = n - 1; j >= i + 1; j--)
			C[i] *= j;
		for (INT32S j = n - 1 - i; j >= 2; j--)
			C[i] /= j;
	}
	for (INT32S k = 0; k < N; k++) {
		FP64 u = (FP64)k / N;
		INT32S x = 0, y = 0;
		for (INT32S i = 0; i < n; i++) {
			FP64 bezBlendFcn = C[i] * pow(u, i) * pow(1 - u, n - 1 - i);
			x += xCtrl[i] * bezBlendFcn;
			y += yCtrl[i] * bezBlendFcn;
		}
		drawPoint(x, y);
	}
}
/* ---------------- DRAW COPY ---------------- */
void Graphics::drawCopy(INT32S x0, INT32S y0, RGBBASIC* gt, INT32S width, INT32S height)
{
	for (INT32S i = 0; i < height; i++) {
		for (INT32S j = 0; j < width; j++) {
			RGB t = *(gt + i * width + j);
			if (t == TRANSPARENT)continue;					//RGB 0xFF000000即透明
			setPoint(x0 + j, y0 + i, t);
		}
	}
}
/* ---------------- FILL ---------------- */
void Graphics::fill(INT32S sx, INT32S sy, INT32S ex, INT32S ey, RGB color)
{
	if (sy > ey) {
		INT32S t = ey;
		ey = sy, sy = t;
	}
	if (sx > ex) {
		INT32S t = ex;
		ex = sx, sx = t;
	}
	for (INT32S y = sy; y <= ey; y++) {
		for (INT32S x = sx; x <= ex; x++) {
			setPoint(x, y, color);
		}
	}
}
/* ---------------- FLOOD FILL ----------------
*	广度优先搜索	队列
** ----------------------------------------*/
void Graphics::fillflood(INT32S x0, INT32S y0, RGB color)
{
	RGB color0 = readPoint(x0, y0);
	INT32S x_step[] = { 0,0,1,-1,1,1,-1,-1 };
	INT32S y_step[] = { 1,-1,0,0,1,-1,1,-1 };
	std::queue<INT32S> Qx, Qy;
	Qx.push(x0);Qy.push(y0);
	setPoint(x0, y0, color);

	while (!Qx.empty()) {
		INT32S x = Qx.front(), y = Qy.front();
		for (INT32S i = 0; i < 8; i++) {
			INT32S xt = x + x_step[i], yt = y + y_step[i];
			if (readPoint(xt, yt) == color0 && !judgeOutRange(xt, yt)) {
				setPoint(xt, yt, color);
				Qx.push(xt); Qy.push(yt);
			}
		}
		Qx.pop(); Qy.pop();
	}
}
/* ---------------- fillTriangle ----------------
*	扫描线填充算法
*	
*	基本思想：
*	每条水平扫描线与多边形的边产生一系列交点，
*	交点之间形成一条一条的线段，该线段上的像素就是需要被填充的像素。
*	将这些交点按照x坐标排序，将排序后的交点两两成对，作为线段的两个端点。
*	水平扫描线从上到下（或从下到上）扫描由多条首尾相连的线段，
*	使用要求的颜色填充该水平线段上的像素。
*
*	利用两张链表列出边，
*	*边表ET: 待相交边的存储表	*活动边表AET: 目前扫描线相交的边表，
*	利用ymin确定什么时候考虑该边，利用ymax确定什么时候放弃该边，利用+dx确定交点.
*
*	算法流程:
*		1.创建链表(链表组): Active-Edge Table活动边表, Edge Table边表组
*		2.计算扫描区域y最大最小值
*		3.基于输入边数据, 初始化ET[y]边表组
*		4.扫描循环开始, 扫描线由Ymin -> Ymax扫描
*			5.将当前ET[y]中边数据，即刚开始相交的边数据插入AET
*			6.AET中边两两配对,填充该扫描线
*			7.删除AET中不再相交的边，即ymax < y+1的边
*			8.更新AET各边同扫描线的相交点x值
** ---------------------------------------- */
struct fillPolygon_Edge{								//边表(链表)
	int ymax;											//ymax:边的下端点
	double x,dx;										//x:当前水平扫描线的交点//dx:斜率m的倒数
	fillPolygon_Edge* next = NULL;
};
void Graphics::fillPolygon(INT32S x[], INT32S y[], INT32S n)
{
	const int ETSzie = 1024;
	fillPolygon_Edge* AET = new fillPolygon_Edge(), *ET[ETSzie];//Active-Edge Table:活动边表//Edge Table边表
	/*------ 计算y最大最小值 ------*/
	int maxY = 0, minY = 0x7FFFFFFF;
	for (int i = 0; i < n; i++) {
		maxY = maxY >= y[i] ? maxY : y[i];
		minY = minY <= y[i] ? minY : y[i];
	}
	for (int i = 0; i <= maxY - minY; i++){
		ET[i] = new fillPolygon_Edge();
	}
	/*------ 建立边表ETy坐标 ------*/
	for (int i = 0; i < n; i++){
		int x1 = x[i], x2 = x[(i + 1) % n];
		int y1 = y[i], y2 = y[(i + 1) % n];
		if (y1 == y2)continue;							//水平线舍弃
		int ymin = y1 < y2 ? y1 : y2;
		fillPolygon_Edge* tE = new fillPolygon_Edge;	//创建新边表节点
		tE->ymax = y1 > y2 ? y1 : y2; 					//下端点ymin,上端点ymax,下端点x,斜率倒数
		tE->x = y1 < y2 ? x1 : x2;
		tE->dx = (double)(x2 - x1) / (y2 - y1); 
		tE->next = ET[ymin - minY]->next;
		ET[ymin - minY]->next = tE;						//插入ET
	}
	/*------ 扫描线由Ymin -> Ymax扫描 ------*/
	for (int y = minY; y <= maxY; y++) {
		fillPolygon_Edge* p = AET;
		/*------ 将当前ET[y]中边数据，即刚开始相交的边数据插入AET ------*/
		while (ET[y - minY]->next) {					//ET[y]链表遍历
			fillPolygon_Edge* pET = ET[y - minY]->next;
			while (p->next) {							//AET链表遍历,在AET中搜索合适的插入位置
				if ((pET->x < p->next->x) || (pET->x == p->next->x && pET->dx < p->next->dx))//按x增序(相等则dx增序)插入AET
					break;								//找到位置
				p = p->next;						
			}
			ET[y - minY]->next = pET->next;				//ET链表中删除pET
			pET->next = p->next; p->next = pET;			//pET插入AET的当前位置
		}
		/*------ AET中边两两配对,填充该扫描线 ------*/
		p = AET;
		while (p->next && p->next->next) {				//链表遍历
			for (int x = p->next->x; x <= p->next->next->x; x++) {
				setPoint(x, y, PaintColor);
			}
			p = p->next->next;
		}
		/*------ 删除AET中不再相交的边 ------*/
		p = AET;
		while (p->next) {
			if (p->next->ymax == y + 1) {
				fillPolygon_Edge* pt = p->next;
				p->next = pt->next;
				delete pt;
			}
			else p = p->next;
		}
		/*------ 更新AET各边同扫描线的相交点x值 ------*/
		p = AET;
		while (p->next) {
			p->next->x += p->next->dx;
			p = p->next;
		}
	}
}
/* ---------------- DRAW CHARACTER ----------------
*	字库: font.h
** ---------------------------------------- */
void Graphics::drawChar(INT32S x0, INT32S y0, CHAR charac)
{
	const INT32S FontLibSize = 16;
	Graphics Maptemp;
	Maptemp.setSize(FontSize, FontSize);
	Maptemp.init();
	Maptemp.clear(TRANSPARENT);
	Maptemp.PaintColor = PaintColor;
	INT32S x = 0, y = 0;
	for (INT32S i = 0; i < FontLibSize; i++) {
		INT8U t = asc2_1608[charac - 32][i];	//调用2412字体
		for (INT32S j = 0; j < 8; j++) {
			if (t & 0x80)Maptemp.drawPoint(x, y);
			t <<= 1;
			y++;
			if (y == FontLibSize) {
				y = 0;	x++;
				break;
			}
		}
	}
	Maptemp.scaling((double)FontSize / FontLibSize, (double)FontSize / FontLibSize);
	if ((double)FontSize / FontLibSize < 1)Maptemp.confirmTrans();
	drawCopy(x0, y0, Maptemp.Map, FontSize, FontSize);
}
/* ---------------- DRAW STRING ---------------- */
void Graphics::drawString(INT32S x0, INT32S y0, const CHAR* str, INT32S n)
{
	for (int i = 0; i < n; i++) {
		drawChar(x0 + FontSize * i, y0, str[i]);
	}
}
/* ---------------- DRAW NUMBER ---------------- */
void Graphics::drawNum(INT32S x0, INT32S y0, FP64 num)
{
	CHAR numstr[100];
	INT32S cur = 0;
	INT32S integer = num;							//整数部分
	FP64 decimal = num - integer;					//小数部分
	if (num < 0) {
		integer *= -1; decimal *= -1;
	}
	/*------ 整数部分 ------*/
	if (integer == 0) numstr[cur++] = '0';
	while (integer > 0) {
		numstr[cur++] = integer % 10 + '0';
		integer /= 10;
	}
	if (num < 0)numstr[cur++] = '-';
	//反转一下
	for (INT32S i = 0; i <= cur / 2; i++) {
		CHAR t = numstr[i];
		numstr[i] = numstr[cur - i - 1]; numstr[cur - i - 1] = t;
	}
	/*------ 小数部分 ------*/
	FP64 accur = 1E-10;								//精度
	if (decimal > accur) {
		numstr[cur++] = '.';
		while (decimal > accur) {
			decimal *= 10;
			numstr[cur++] = (CHAR)decimal + '0';
			decimal -= (INT32S)decimal;
		}
	}
	drawString(x0, y0, numstr, cur);
}
/******************************************************************************
*                    二维变换
******************************************************************************/
/* ---------------- TRANSLATION ---------------- */
void Graphics::translation(INT32S dx, INT32S dy)
{
	Mat<FP64> M(3);
	M(0, 2) = dx; M(1, 2) = dy;
	gM.mult(M, gM);
}
/* ---------------- ROMOTE ---------------- */
void Graphics::rotate(FP64 theta, INT32S x0, INT32S y0) 
{
	translation(-x0, -y0);
	Mat<FP64> M(3);
	M(0, 0) = cos(theta); M(0, 1) = -1 * sin(theta);
	M(1, 0) = sin(theta); M(1, 1) = cos(theta);
	gM.mult(M, gM);
	translation(x0, y0);
}
/* ---------------- scaling ---------------- */
void Graphics::scaling(FP64 sx, FP64 sy)
{
	Mat<FP64> M(3);
	M(0, 0) = sx; M(1, 1) = sy;
	if (sx <= 1 && sy <= 1) { gM.mult(M, gM); return; }
	Graphics Maptemp;
	Maptemp.setSize(gWidth, gHeight);
	Maptemp.init();
	Maptemp.clear(TRANSPARENT);
	Mat<FP64> p0(3, 1), p1(3, 1), p2(3, 1);
	M(2, 0) = 1;
	for (INT32S y = 0; y < gHeight; y++) {
		for (INT32S x = 0; x < gWidth; x++) {
			p0(0, 0) = (FP64)x;	p0(1, 0) = (FP64)y;
			p1.mult(M, p0);
			RGB t = readPoint(x, y);
			if (t == TRANSPARENT)continue;
			p0(0, 0) = (FP64)x + 1;	p0(1, 0) = (FP64)y + 1;
			p2.mult(M, p0);
			Maptemp.fill(p1(0, 0), p1(1, 0), p2(0, 0), p2(1, 0), t);
		}
	}
	free(Map);
	Map = Maptemp.Map;	Maptemp.Map = NULL;
}
/******************************************************************************
*                    SET 设置
******************************************************************************/
/* ---------------- setSize ---------------- */
void Graphics::setSize(INT32S width, INT32S height)
{
	gWidth = width;	gHeight = height;
}
/* ---------------- judgeOutRange ---------------- */
bool Graphics::judgeOutRange(INT32S x0, INT32S y0)
{
	if (x0<0 || x0>= gWidth)return true;
	if (y0<0 || y0>= gHeight)return true;
	return false;
}