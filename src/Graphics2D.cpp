#include "Graphics2D.h"

/********************************************************************
*
*                    2维 图形学
*
********************************************************************/

ARGB Graphics::PaintColor = 0xFFFFFFFF;
int  Graphics::PaintSize = 0, Graphics::FontSize = 16;

/*---------------- 画点 ----------------*/
void Graphics::drawPoint(Mat<ARGB>& image, int x0, int y0)
{
	if (image.isOut(x0, y0)) return;
	image(x0, y0) = PaintColor;										//基础点(点粗==0)

	// 点粗>0时
	if (PaintSize > 0) {
		int x = 0, y = PaintSize, p = 3 - (PaintSize << 1);			//初始点:天顶(0,r)//p:决策参数(r右移即乘2)
		int x_step[] = { 1,1,-1,-1 },
			y_step[] = { 1,-1,1,-1 };								//上下左右对称四个点

		// 绘制圆 (x=0始,y=x终)
		while (x <= y) {
			for (int i = x0 - x; i <= x0 + x; i++) { 				//填充圆内
				if (!image.isOut(i, y0 - y))
					image(i, y0 - y) = PaintColor;
				if (!image.isOut(i, y0 + y))
					image(i, y0 + y) = PaintColor;
			}

			for (int i = x0 - y; i <= x0 + y; i++) { 
				if (!image.isOut(i, y0 - x))
					image(i, y0 - x) = PaintColor;
				if (!image.isOut(i, y0 + x))
					image(i, y0 + x) = PaintColor;
			}

			x++;
			int dp = 4 * x + 6;
			if (p < 0) p += dp;
			else { p += dp - 4 * y + 4; y--; }							//p过临界值，y加一
		}
	}
}

/*---------------- 画线 : Bresenham 算法 ----------------*/
void Graphics::drawLine(Mat<ARGB>& image, int x1, int y1, int x2, int y2)
{
	int inc[2] = { 0 },
		delta[2] = { x2 - x1, y2 - y1 },
		p[2] = { x1, y1 };

	for (int d = 0; d < 2; d++) {
		inc[d] = delta[d] == 0 ? 0 : (delta[d] > 0 ? 1 : -1);
		delta[d] *= inc[d];    // |Δ| 
	}

	int a = 0, d1, d2;

	if (delta[0] >= delta[1]) {
		d1 = 0;
		d2 = 1;
	}
	else {
		d1 = 1;
		d2 = 0;
	}

	for (int i = 0; i <= delta[d1]; i++) {
		drawPoint(image, p[0], p[1]);

		p[d1] += inc[d1];

		a += delta[d2];
		if (a >= delta[d1]) {
			a -= delta[d1];
			p[d2] += inc[d2];
		}
	}
}

/*---------------- 画折线 ----------------*/
void Graphics::drawLine(Mat<ARGB>& image, int* x, int* y, int n) {
	int xt, yt;
	for (int i = 0; i < n; i++) {
		if (i != 0)
			drawLine(image, xt, yt, x[i], y[i]);
		xt = x[i];
		yt = y[i];
	}
}

/*---------------- 画三角 ----------------*/
void Graphics::drawTriangle(Mat<ARGB>& image, 
							int x1, int y1, int x2, int y2, int x3, int y3)
{
	drawLine(image, x1, y1, x2, y2);
	drawLine(image, x2, y2, x3, y3);
	drawLine(image, x3, y3, x1, y1);
}

/*---------------- 画矩形 ----------------*/
void Graphics::drawRectangle(Mat<ARGB>& image, int x1, int y1, int x2, int y2)
{
	drawLine(image, x1, y1, x2, y1);
	drawLine(image, x1, y1, x1, y2);
	drawLine(image, x1, y2, x2, y2);
	drawLine(image, x2, y1, x2, y2);
}

void Graphics::drawRegularPolygon (Mat<ARGB>& image, int  x, int  y, int l, int n, double a0) {
    double a = 2 * PI / n;

    for (int i = 1; i < n; i++) 
		drawLine(image, 
			x + l * cos(a0 + a * (i-1)), y + l * sin(a0 + a * (i-1)),
			x + l * cos(a0 + a * i),	 y + l * sin(a0 + a * i));

	drawLine(image, 
		x + l * cos(a0 + a * (n-1)), y + l * sin(a0 + a * (n-1)),
		x + l * cos(a0), y + l * sin(a0));
}

void Graphics::drawPolygon(Mat<ARGB>& image, int* x, int* y, int n)
{
	for (int i = 0; i < n; i++)
		drawLine(image, x[i], y[i], x[(i + 1) % n], y[(i + 1) % n]);
}



/*---------------- 画圆 : Bresenham 算法 ---------------*/
void Graphics::drawCircle(Mat<ARGB>& image, int x0, int y0, int r)
{
	int x = 0, y = r, p = 3 - (r << 1);  //初始点:天顶(0,r)
	int x_step[] = { 1, 1,-1,-1 }, 
		y_step[] = { 1,-1, 1,-1 };		//上下左右对称四个点

	while (x <= y) {
		for (int i = 0; i < 4; i++) {  //八分圆对称绘制
			drawPoint(image, x0 + x * x_step[i], y0 + y * y_step[i]);
			drawPoint(image, x0 + y * y_step[i], y0 + x * x_step[i]);
		}

		x++;
		int dp = 4 * x + 6;

		if (p < 0)
			p += dp;
		else { 
			p += dp - 4 * y + 4; 
			y--; 
		} 
	}
}

/*---------------- 画椭圆 ---------------*/
void Graphics::drawEllipse(Mat<ARGB>& image, int x0, int y0, int rx, int ry)
{
	int64 rx2 = rx * rx, ry2 = ry * ry;
	int x = 0, y = ry;											//初始点:天顶
	int64 p = ry2 + rx2 * (0.25 - ry);								//p:决策参数
	int x_step[] = { 1,1,-1,-1 }, y_step[] = { 1,-1,1,-1 };		//上下左右对称四个点

	// 绘制椭圆 (区域1)
	while (ry2 * x < rx2 * y) {
		for (int i = 0; i < 4; i++) 								//四分对称绘制
			drawPoint(image, x0 + x * x_step[i], y0 + y * y_step[i]);
		x++;
		int64 dp = (1 + 2 * x) * ry2;
		if (p < 0)p += dp;
		else { y--; p += dp - 2 * rx2 * y; }						//p过临界值，y加一
	}

	// 绘制椭圆 (区域2)
	p = ry2 * (x + 0.5) * (x + 0.5) + rx2 * (y - 1) * (y - 1) - ry2 * rx2;
	while (y >= 0) {
		for (int i = 0; i < 4; i++) 								//八分圆对称绘制
			drawPoint(image, x0 + x * x_step[i], y0 + y * y_step[i]);
		y--;
		int64 dp = (1 - 2 * y) * rx2;
		if (p > 0)p += dp;
		else { x++; p += dp + 2 * ry2 * x; }						//p过临界值，y加一
	}
}

/*---------------- 画网格 ----------------*/
void Graphics::drawGrid(Mat<ARGB>& image, int sx, int sy, int ex, int ey, int dx, int dy)
{
	if (dx > 0) {
		for (int x = sx; x <= ex; x += dx) {
			if (image.isOut(x, 0)) continue;
			drawLine(image, x, sy, x, ey);
		}
	}
	else {
		for (int x = sx; x >= ex; x += dx) {
			if (image.isOut(x, 0)) continue;
			drawLine(image, x, sy, x, ey);
		}
	}
	if (dy > 0) {
		for (int y = sy; y <= ey; y += dy) {
			if (image.isOut(0, y)) continue;
			drawLine(image, sx, y, ex, y);
		}
	}
	else {
		for (int y = sy; y >= ey; y += dy) {
			if (image.isOut(0, y)) continue;
			drawLine(image, sx, y, ex, y);
		}
	}
}

/*----------------[ DRAW BEZIER CURVE ]----------------*/
void Graphics::drawBezier(Mat<ARGB>& image, vector<vector<double>>& points, int n)
{
	vector<double> p;

	for (int i = 0; i <= n; i++) {
		double t = i / (double)n;
		BezierCurve(points, n, p);
		drawPoint(image, (int)p[0], (int)p[1]);
	}
}

/*---------------- 填充 ----------------*/
void Graphics::fillRectangle(Mat<ARGB>& image, int sx, int sy, int ex, int ey, ARGB color)
{
	if (sy > ey) std::swap(ey, sy);
	if (sx > ex) std::swap(ex, sx);

	for (int y = sy; y <= ey; y++)
		for (int x = sx; x <= ex; x++)
			image(x, y) = color;
}

/*---------------- FLOOD 填充 ----------------*/
void Graphics::fillFlood(Mat<ARGB>& image, int x0, int y0, ARGB color)
{
	ARGB color0 = image(x0, y0);
	int x_step[] = { 0,0,1,-1,1,1,-1,-1 },
		y_step[] = { 1,-1,0,0,1,-1,1,-1 };
	std::queue<int> Qx, Qy;
	Qx.push(x0);
	Qy.push(y0);
	image(x0, y0) = color;

	while (!Qx.empty()) {
		int x = Qx.front(),
			y = Qy.front();
		for (int i = 0; i < 8; i++) {
			int xt = x + x_step[i],
				yt = y + y_step[i];
			if (image(xt, yt) == color0 && !image.isOut(xt, yt)) {
				image(xt, yt) = color;
				Qx.push(xt);
				Qy.push(yt);
			}
		}
		Qx.pop();
		Qy.pop();
	}
}

/*---------------- 扫描线填充 ----------------*/
struct fillPolygon_Edge {								//边表(链表)
	int ymax;											//ymax:边的下端点
	double x, dx;										//x:当前水平扫描线的交点//dx:斜率m的倒数
	fillPolygon_Edge* next = NULL;
};

void Graphics::fillPolygon(Mat<ARGB>& image, int* x, int* y, int n)
{
	const int ETSzie = 1024;
	fillPolygon_Edge* AET = new fillPolygon_Edge(), * ET[ETSzie];//Active-Edge Table:活动边表//Edge Table边表

	// 计算y最大最小值 ------
	int maxY = 0, minY = 0x7FFFFFFF;
	for (int i = 0; i < n; i++) {
		maxY = maxY >= y[i] ? maxY : y[i];
		minY = minY <= y[i] ? minY : y[i];
	}
	for (int i = 0; i <= maxY - minY; i++)
		ET[i] = new fillPolygon_Edge();

	// 建立边表ETy坐标 ------
	for (int i = 0; i < n; i++) {
		int x1 = x[i], x2 = x[(i + 1) % n],
			y1 = y[i], y2 = y[(i + 1) % n];

		if (y1 == y2) continue;							//水平线舍弃

		int ymin = y1 < y2 ? y1 : y2;

		fillPolygon_Edge* tE = new fillPolygon_Edge;	//创建新边表节点
		tE->ymax = y1 > y2 ? y1 : y2; 					//下端点ymin,上端点ymax,下端点x,斜率倒数
		tE->x = y1 < y2 ? x1 : x2;
		tE->dx = (double)(x2 - x1) / (y2 - y1);
		tE->next = ET[ymin - minY]->next;
		ET[ymin - minY]->next = tE;						//插入ET
	}

	// 扫描线由Ymin -> Ymax扫描
	for (int y = minY; y <= maxY; y++) {
		fillPolygon_Edge* p = AET;
		// 将当前ET[y]中边数据，即刚开始相交的边数据插入AET
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
		// AET中边两两配对,填充该扫描线
		p = AET;
		while (p->next && p->next->next) {				//链表遍历
			for (int x = p->next->x; x <= p->next->next->x; x++)
				if(!image.isOut(x, y))
					image(x, y) = PaintColor;
			p = p->next->next;
		}
		// 删除AET中不再相交的边
		p = AET;
		while (p->next) {
			if (p->next->ymax == y + 1) {
				fillPolygon_Edge* pt = p->next;
				p->next = pt->next;
				delete pt;
			}
			else p = p->next;
		}
		// 更新AET各边同扫描线的相交点x值
		p = AET;
		while (p->next) {
			p->next->x += p->next->dx;
			p = p->next;
		}
	}
}

/*---------------- 画字符 ----------------*/
void Graphics::drawChar(Mat<ARGB>& image, int x0, int y0, char charac)
{
	static const int FontLibSize = 16;
	int k = FontSize / FontLibSize + 1;
	int x = 0, y = 0;

	for (int i = 0; i < FontLibSize; i++) {
		INT8U t = asc2_1608[charac - 32][i];

		for (int j = 0; j < 8; j++) {
			if (t & 0x80)
				fillRectangle(
					image,
					x * k + x0, y * k + y0, 
					(x + 1) * k + x0, 
					(y + 1) * k + y0, 
					PaintColor
				);
			t <<= 1; 
			x++;

			if (x == FontLibSize) { 
				x = 0;	
				y++; 
				break; 
			}
		}
	}
}
/*---------------- 画字符串 ----------------*/
void Graphics::drawString(Mat<ARGB>& image, int x0, int y0, const char* str)
{
	for (int i = 0; str[i] != '\0'; i++) 
		drawChar(image, x0, y0 + FontSize * i, str[i]);
}

/*---------------- 画数字 ----------------*/
void Graphics::drawNum(Mat<ARGB>& image, int x0, int y0, fp64 num)
{
	char numstr[100];
	int cur = 0;
	int integer = num;							//整数部分
	fp64 decimal = num - integer;					//小数部分
	if (num < 0)
		integer *= -1, decimal *= -1;

	// 整数部分
	if (integer == 0) 
		numstr[cur++] = '0';

	while (integer > 0) {
		numstr[cur++] = integer % 10 + '0';
		integer /= 10;
	}

	if (num < 0) 
		numstr[cur++] = '-';

	for (int i = 0; i < cur / 2; i++)			//反转一下
		std::swap(numstr[i], numstr[cur - i - 1]);

	// 小数部分
	fp64 accur = 1E-10;								//精度
	if (decimal > accur) {
		numstr[cur++] = '.';
		while (decimal > accur) {
			decimal *= 10;
			numstr[cur++] = (char)decimal + '0';
			decimal -= (int)decimal;
		}
	} 
	numstr[cur] = '\0';
	drawString(image, x0, y0, numstr);
}
