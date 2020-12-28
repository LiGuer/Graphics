#include "LiGu_Graphics/Graphics.h"

int main()
{
	// init
	Graphics g(2000, 2000);
	// Draw Point
	g.PaintColor = 0xFFFFFF;
	for (int i = 0; i < 50; i++) { g.PaintSize = i; g.drawPoint(50 + i * 38, 1800); }
	// Draw Circle & Ellipse
	g.PaintSize = 5; g.PaintColor = 0x66CCFF; g.drawEllipse(1000, 1000, 500, 300);
	g.PaintColor = 0xFF0000; g.drawCircle(1000, 1000, 500);
	g.PaintSize = 3; g.PaintColor = 0xCCFF66; g.drawEllipse(1000, 1000, 200, 500);
	// Draw & Fill drawRectangle
	g.PaintColor = 0xCCFF00; g.drawRectangle(1900, 100, 1700, 300);
	g.fillRectangle(1850, 150, 1750, 250, 0xCCFF00);
	{ // Draw Polygon
		int xt[] = { 100 ,300,100 }, yt[] = { 900 ,1500,1700 };
		g.PaintColor = 0xCC0066; g.drawPolygon(xt, yt, 3);
	}
	{ // Draw Wave
		const int N = 100;
		int x[N], y[N];
		for (int i = 0; i < N; i++) {
			x[i] = i * 10 + 500; y[i] = sin(0.02 * x[i]) * 100 + 200;
		}
		g.PaintColor = 0xFF0000; g.drawWave(x, y, N);
	}
	// Draw Num String & String
	g.drawNum(1000, 1200, -1203.567);
	g.PaintColor = 0x0099FF;  g.FontSize = 100; g.drawString(0, 0, "Ligu~", 5);
	{ // Draw Bezier Curve and Fill by Flood
		g.PaintColor = 0x00EEFF; g.PaintSize = 3;
		int x[] = { 1900 - 200,1900,1900 + 200,1900 ,1900 - 200 };
		int y[] = { 600,600 - 500,600,600 - 100 ,600 };
		g.drawBezier(x, y, 5);
		g.fillFlood(1900, 502, 0x00EEFF);
	}
	{ // Love Heart
		g.PaintColor = 0xFF5500; g.PaintSize = 3;
		int x[] = { 1800,1800 + 100,1800 + 250 ,1800 };
		int y[] = { 1000,1000 - 150,1000  ,1000 + 230 };
		g.drawBezier(x, y, 4);
		int x2[] = { 1800 ,1800 - 100 ,1800 - 250  ,1800 };
		int y2[] = { 1000,1000 - 150,1000  ,1000 + 230 };
		g.drawBezier(x2, y2, 4);
		g.FontSize = 30; g.drawString(1700, 1150, "I Love U", 8);
	}
	{ // Rotate & Scaling & Translation
		g.rotate(3.14 * 60 / 360, 500, 500);
		g.scaling(0.2, 0.2, 0, 0);
		g.translation(100, 50);
		g.PaintSize = 5; g.PaintColor = 0x66CCFF; g.drawEllipse(1000, 1000, 500, 300);
		g.PaintColor = 0xFF0000; g.drawCircle(1000, 1000, 500);
		g.PaintSize = 3; g.PaintColor = 0xCCFF66; g.drawEllipse(1000, 1000, 200, 500);
		g.PaintColor = 0xCCFF00; g.drawRectangle(1900, 100, 1600, 400);
		g.fillRectangle(1800, 200, 1700, 300, 0xCCFF00);
		g.TransMat.E(3);
	}
	{ // Translucent
		int sx = 280, sy = 250;
		g.fillRectangle(sx, sy, sx + 200, sy + 200, 0x8800FFFF);
		g.fillRectangle(sx + 500, sy + 300, sx + 700, sy + 500, 0x88FF00FF);
		g.fillRectangle(sx + 800, sy + 800, sx + 1100, sy + 1100, 0x88FFFF00);
		g.fillRectangle(sx + 800, sy + 1400, sx + 1100, sy + 1700, 0x8800FF00);
		g.fillRectangle(sx + 100, sy + 800, sx + 400, sy + 1100, 0x88FF0000);
		g.fillRectangle(sx + 500, sy + 800, sx + 700, sy + 1100, 0x880000FF);
	}
	{ // Sub Graphics
		Graphics gt(g.gWidth, g.gHeight);
		gt.drawCopy(0, 0, g.Map, g.gWidth, g.gHeight);
		gt.scaling(0.2, 0.2, 0, 0);
		gt.transSelf();
		gt.CutSelf(0, 0, g.gWidth / 5, g.gHeight / 5);
		g.drawCopy(0, 500, gt.Map, gt.gWidth, gt.gHeight);
	}
	// Save .PPM
	g.PicWrite("D:/LIGU.ppm");
}