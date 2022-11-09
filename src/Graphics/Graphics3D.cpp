#include "Graphics3D.h"

extern vector<double> Graphics::lightVector = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3) };

bool Graphics::drawPoint (Mat<ARGB>& image, Mat<int>& Z_buf, int x, int y, int z) 
{
    if (image.isOut(x, y) || z <= Z_buf(x, y))
        return false;

    image(x, y) = PaintColor;
	Z_buf(x, y) = z; 

    return true;
}

bool Graphics::drawPoint(Mat<ARGB>& image, Mat<int>& Z_buf, int x, int y, int z, double fx, double fy, double fz)
{
    if (image.isOut(x, y) || z <= Z_buf(x, y))
        return false;

    double a = (fx * lightVector[0] + fy * lightVector[1] + fz * lightVector[2]) / sqrt(fx * fx + fy * fy + fz * fz);

    if (a < 0)
        a = 0;

    image(x, y) = ((ARGB)(a * (unsigned char) PaintColor)) +
                  ((ARGB)(a * (unsigned char)(PaintColor >> 8 )) << 8 ) +
                  ((ARGB)(a * (unsigned char)(PaintColor >> 16)) << 16);
    Z_buf(x, y) = z;

    return true;
}

void Graphics::drawLine (Mat<ARGB>& image, Mat<int>& Z_buf, 
    int sx, int ex, 
    int sy, int ey, 
    int sz, int ez
) {
	int a    [3] = { 0 },
		inc  [3] = { 0 },
		delta[3] = { ex - sx, ey - sy, ez - sz },
		p    [3] = { sx, sy, sz };
 
	for (int d = 0; d < 3; d++) {
		inc[d] = delta[d] == 0 ? 0 : (delta[d] > 0 ? 1 : -1);
		delta[d] *= inc[d];    // |Δ| 
	}
	int max_delta = max(max(delta[0], delta[1]), delta[2]);

	for (int i = 0; i <= max_delta; i++) {
		drawPoint(image, Z_buf, p[0], p[1], p[2]);

		for (int d = 0; d < 3; d++) {
			a[d] += delta[d];
			if (a[d] >= max_delta) {
				a[d] -= max_delta;
				p[d] += inc[d];
			}
		}
	}
}


void Graphics::drawLine(Mat<ARGB>& image, Mat<int>& Z_buf, vector<vector<int>>& p, bool close) {
    int n = p.size();

	for (int i = 0; i < n - 1; i++) 
        drawLine(image, Z_buf, 
            p[i][0], p[i + 1][0], 
            p[i][1], p[i + 1][1], 
            p[i][2], p[i + 1][2]
        );

	if (close) 
        drawLine(image, Z_buf, 
            p[0][0], p[n - 1][0],
            p[0][1], p[n - 1][1],
            p[0][2], p[n - 1][2]
        );
}

/*
void Graphics::drawTriangle(Mat<ARGB>& image, Mat<int>& Z_buf, vector<int>& p1, vector<int>& p2, vector<int>& p3) 
{

    int pXminI = 0;

    for (int i = 1; i < 3; i++) 
        pXminI = pt[i][0] < pt[pXminI][0] ? i : pXminI;

    std::swap(pt[0], pt[pXminI]);

    if(pt[0][0] >= g.Canvas.rows)
        return;
    if(std::max(pt[1][0], pt[2][0]) < 0)
        return;
    if(std::min(pt[0][1], std::min(pt[1][1], pt[2][1])) >=g.Canvas.cols)
        return;
    if(std::max(pt[0][1], std::max(pt[1][1], pt[2][1])) < 0)
        return;

    //[3]
    static Mat<int> err[2], inc[2], delta[2], point[2];

    for (int k = 0; k < 2; k++) {
        err[k].zero(Dim);
        inc[k].zero(Dim);
        delta[k].sub(pt[k + 1], pt[0]); 
        point[k] = pt[0];
        for (int dim = 0; dim < Dim; dim++) {
            inc  [k][dim] = delta[k][dim] == 0 ? 0 : (delta[k][dim] > 0 ? 1 : -1);//符号函数(向右,垂直,向左)
            delta[k][dim] = abs(delta[k][dim]);	
        }
    }

    int dXMaxCur = delta[0][0] > delta[1][0] ? 0 : 1;

    bool flag = true;															//三角形转折点检测开关(不然会二次检测)

    for (int i = 0; i <= delta[dXMaxCur][0]; i++) {

        //[4]三角形转折点检测
        if (i == delta[1 - dXMaxCur][0] && flag) {
            flag = false;														//关闭检测
            int kt = 1 - dXMaxCur;
            delta[kt].sub(pt[dXMaxCur + 1], pt[kt + 1]);
            point[kt] = pt[kt + 1];
            for (int dim = 1; dim < Dim; dim++) {
                inc  [kt][dim]  = delta[kt][dim] == 0 ? 0 : (delta[kt][dim] > 0 ? 1 : -1);	//符号函数(向右,垂直,向左)
                delta[kt][dim] *= delta[kt][dim] < 0 ? -1 : 1;					//向左
            }
        }

        //[5]画线
        static Mat<int> errTmp(Dim), incTmp(Dim), deltaTmp, pointTmp; 
        errTmp.zero(); 
        pointTmp = point[0];

        deltaTmp.sub(point[1], point[0]);

        for (int dim = 0; dim < Dim; dim++) {									//设置xyz单步方向	
            incTmp  [dim] = deltaTmp[dim] == 0 ? 0 : (deltaTmp[dim] > 0 ? 1 : -1);//符号函数(向右,垂直,向左)
            deltaTmp[dim] = abs(deltaTmp[dim]);	
        }

        int distanceTmp = deltaTmp.max();										//总步数

        for (int i = 0; i <= distanceTmp; i++) {								//画线
            setPix(pointTmp, 0, FaceColorTmp); 

            for (int dim = 0; dim < Dim; dim++) {								//xyz走一步
                errTmp[dim] += deltaTmp[dim];
                if (errTmp  [dim] >= distanceTmp) { 
                    errTmp  [dim] -= distanceTmp; 
                    pointTmp[dim] += incTmp[dim];
                }
            }
        }

        for (int k = 0; k < 2; k++) {
            point[k][0]++;
            for (int dim = 1; dim < Dim; dim++) {								//xyz走一步
                err[k][dim] += delta[k][dim];
                if (delta[k][0] == 0) break;
                if (err[k][dim] >= delta[k][0]) {
                    point[k][dim] += err[k][dim] / delta[k][0] * inc[k][dim]; 
                    err  [k][dim]  = err[k][dim] % delta[k][0];
                }
            }
        }

    }
}

void Graphics::drawTriangleSet(Mat<ARGB>& image, Mat<int>& Z_buf, vector<vector<vector<int>>>& p) {
    int n = p.size();

	for (int i = 0; i < n; i++) 
		drawTriangle(
            image, Z_buf,
			p[i][0], 
			p[i][1], 
			p[i][2]
		);
}*/


void Graphics::drawSphere(Mat<ARGB>& image, Mat<int>& Z_buf, int x0, int y0, int z0, int r)
{
    int x_step[] = { 1, 0, 1, 1, 0, 1 },
        y_step[] = { 0, 1, 1, 0, 1, 1 },
        z_step[] = { 0, 0, 0,-1,-1,-1 },
        sign[3][8] = { 
            { 1,-1, 1, 1,-1, 1,-1,-1 },
            { 1, 1,-1, 1,-1,-1, 1,-1 },
            { 1, 1, 1,-1, 1,-1,-1,-1 } 
        },
        permu[6][3] = {
            {0, 1, 2}, {0, 2, 1}, 
            {1, 0, 2}, {1, 2, 0},
            {2, 0, 1}, {2, 1, 0},
        },
        e[6];
    queue<vector<int>> Q;
    map<vector<int>, int> M;

    Q.push({ 0, 0, r });

    while (!Q.empty()) {
        vector<int> p = Q.front();
        Q.pop();
        {
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 6; j++) {
                    drawPoint(image, Z_buf,
                        x0 + p[permu[j][0]] * sign[permu[j][0]][i],
                        y0 + p[permu[j][1]] * sign[permu[j][1]][i],
                        z0 + p[permu[j][2]] * sign[permu[j][2]][i],
                        p[permu[j][0]] * sign[permu[j][0]][i],
                        p[permu[j][1]] * sign[permu[j][1]][i],
                        p[permu[j][2]] * sign[permu[j][2]][i]
                    );
                }
            }
        }

        for (int i = 0; i < 6; i++) {
            int x = p[0] + x_step[i],
                y = p[1] + y_step[i],
                z = p[2] + z_step[i];

            if (x > y || y > z || M.find({ x, y, z }) != M.end()) {
                e[i] = DBL_MAX;
                continue;
            }

            e[i] = abs(x * x + y * y + z * z - r * r);

            if (e[i] <= r) {
                Q.push({ x, y, z });
                M[{x, y, z}] = 1;
            }
        }

    }

}


void Graphics::drawFunction(Mat<ARGB>& image, Mat<int>& Z_buf, int xs, int ys, int zs, 
    function<double(double, double, double)> f, 
    function<void(double, double, double, double&, double&, double&)> df)
{

    int d[26][3];
    double e, fx, fy, fz;
    queue<vector<int>> Q;
    map<vector<int>, int> M;

    for (int i = 1; i < 27; i++) {
        d[i - 1][0] = (i % 3) - 1;
        d[i - 1][1] =((i % 9) / 3) - 1;
        d[i - 1][2] = (i / 9) - 1;
    }

    Q.push({ xs, ys, zs });
    M[{ xs, ys, zs }] = 1;

    while (!Q.empty()) {
        vector<int> p = Q.front();
        Q.pop();

        {
            df(p[0], p[1], p[2], fx, fy, fz);
            drawPoint(image, Z_buf,
                p[0], p[1], p[2], fx, fy, fz
            );
        }

        for (int i = 0; i < 26; i++) {
            int x = p[0] + d[i][0],
                y = p[1] + d[i][1],
                z = p[2] + d[i][2];

            if (x < 0 || x >= image.rows ||
                y < 0 || y >= image.cols ||
                M.find({ x, y, z }) != M.end()) {
                e = DBL_MAX;
                continue;
            }
            M[{x, y, z}] = 1;

            e= f(x, y, z);

            if (e <= 0) {
                Q.push({ x, y, z });
            }
        }

    }
}