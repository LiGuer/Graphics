#include "Graphics.h"
/******************************************************************************
*                    Basic Function
******************************************************************************/
/* ---------------- INIT ---------------- */
void Graphics::init() {
	if (Map != NULL)free(Map);
	Map = (RGB*)malloc(sizeof(RGB) * gWidth * gHeight);
	clear(0);
	gM.E(3);
}
/* ---------------- CLEAR ---------------- */
void Graphics::clear(RGB color)
{
	if (color == TRANSPARENT || color == 0) {	//memset按字节处理，故只能处理高低字节相同的值
		memset(Map, color, sizeof(RGB) * gWidth * gHeight);
		return;
	}
	for (INT32S y = 0; y < gHeight; y++) {
		for (INT32S x = 0; x < gWidth; x++) {
			setPoint(x, y, color);
		}
	}
}
/* ---------------- SET/READ POINT ---------------- */
void Graphics::setPoint(INT32S x, INT32S y,RGB color) {
	if (judgeOutRange(x, y))return;
	*(Map + y * gWidth + x) = color;
}

RGB Graphics::readPoint(INT32S x, INT32S y) {
	if (judgeOutRange(x, y))return TRANSPARENT;
	return *(Map + y * gWidth + x);
}
/* ---------------- PicWrite ---------------- */
void Graphics::PicWrite(const CHAR* filename) {		// 太慢
	FILE* fp = fopen(filename, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", gWidth, gHeight);// 写图片格式、宽高、最大像素值

	unsigned char color;
	for (INT32S i = 0; i < gHeight; i++) {
		for (INT32S j = 0; j < gWidth; j++) {
			for (INT32S k = 0; k < 3; k++) {
				color = readPoint(j, i) >> (8 * k);
				fwrite(&color, sizeof(color), 1, fp);// 写RGB数据
			}
		}
	}fclose(fp);
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
	p.setValue(2, 0, 1);
	for (INT32S y = 0; y < gHeight; y++) {
		for (INT32S x = 0; x < gWidth; x++) {
			p.setValue(0, 0, (FP64)x); p.setValue(1, 0, (FP64)y);
			p.mult(gM, p, p);
			RGB t = readPoint(x, y);
			if (t != TRANSPARENT)Maptemp.setPoint(p.getValue(0, 0), p.getValue(1, 0), t);
			if (t != TRANSPARENT)Maptemp.setPoint(p.getValue(0, 0) + 0.4, p.getValue(1, 0) + 0.4, t);	//#补丁
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
	for (INT32S r = 1; r <= PaintSize; r++) {
		INT32S x = 0, y = r, p = 3 - (r << 1);		//初始点:天顶(0,r)//p:决策参数(r右移即乘2)
		INT32S x_step[] = { 1,1,-1,-1 }, y_step[] = { 1,-1,1,-1 };		//上下左右对称四个点
		/*------ 绘制圆 (x=0始,y=x终) ------*/
		while (x <= y) {
			for (int i = 0; i < 4; i++) {
				setPoint(x0 + x * x_step[i], y0 + y * y_step[i], PaintColor);
				setPoint(x0 + y * x_step[i], y0 + x * y_step[i], PaintColor);
				setPoint(x0 + x * x_step[i], y0 + (y - 1) * y_step[i], PaintColor);	//#补丁，需优化
				setPoint(x0 + (y - 1) * x_step[i], y0 + x * y_step[i], PaintColor);
			}
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
void Graphics::drawTriangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2, INT32S x3, INT32S y3)
{
	drawLine(x1, y1, x2, y2);
	drawLine(x1, y1, x3, y3);
	drawLine(x2, y2, x3, y3);
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
void Graphics::drawCopy(INT32S x0, INT32S y0, RGB* gt, INT32S width, INT32S height)
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
void Graphics::floodfill(INT32S x0, INT32S y0, RGB color) 
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
/* ---------------- DRAW CHARACTER ---------------- 
*	字库: font.h
** ---------------------------------------- */
void Graphics::drawChar(INT32S x0, INT32S y0, CHAR charac)
{
	INT32S x = x0, y = y0;
	for (INT32S i = 0; i < FontSize; i++){
		INT8U temp = asc2_1608[charac - 32][i];	//调用2412字体
		for (INT32S j = 0; j < 8; j++){
			if (temp & 0x80)drawPoint(x, y);
			temp <<= 1;
			y++;
			if ((y - y0) == FontSize){
				y = y0;	x++;
				break;
			}
		}
	}
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
	Mat<FP64> M;
	M.E(3);
	M.setValue(0, 2, dx);	M.setValue(1, 2, dy);
	gM.mult(M, gM, gM);
}
/* ---------------- ROMOTE ---------------- */
void Graphics::rotate(FP64 theta, INT32S x0, INT32S y0) 
{
	translation(-x0, -y0);
	Mat<FP64> M;
	M.E(3);
	M.setValue(0, 0, cos(theta));	M.setValue(0, 1, -1 * sin(theta));
	M.setValue(1, 0, sin(theta));	M.setValue(1, 1, cos(theta));
	gM.mult(M, gM, gM);
	translation(x0, y0);
}
/* ---------------- scaling ---------------- */
void Graphics::scaling(FP64 sx, FP64 sy)
{
	Mat<FP64> M;
	M.E(3);
	M.setValue(0, 0, sx);	M.setValue(1, 1, sy);
	if (sx <= 1 && sy <= 1) {
		gM.mult(M, gM, gM);
		return;
	}
	Graphics Maptemp;
	Maptemp.setSize(gWidth, gHeight);
	Maptemp.init();
	Maptemp.clear(TRANSPARENT);
	Mat<FP64> p0, p1, p2;
	p0.zero(3, 1); p1.zero(3, 1); p2.zero(3, 1);
	p0.setValue(2, 0, 1);
	for (INT32S y = 0; y < gHeight; y++) {
		for (INT32S x = 0; x < gWidth; x++) {
			p0.setValue(0, 0, (FP64)x); p0.setValue(1, 0, (FP64)y);
			p1.mult(M, p0, p1);
			RGB t = readPoint(x, y);
			if (t == TRANSPARENT)continue;
			p0.setValue(0, 0, (FP64)x + 1); p0.setValue(1, 0, (FP64)y + 1);
			p2.mult(M, p0, p2);
			Maptemp.fill(p1.getValue(0, 0), p1.getValue(1, 0), p2.getValue(0, 0), p2.getValue(1, 0), t);
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
BOOL Graphics::judgeOutRange(INT32S x0, INT32S y0)
{
	if (x0<0 || x0>= gWidth)return true;
	if (y0<0 || y0>= gHeight)return true;
	return false;
}