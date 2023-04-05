#include "Graphics3D.h"

bool Graphics::drawPoint (Mat<ARGB>& image, Mat<int>& Z_buf, int x, int y, int z) 
{
    if (image.isOut(x, y) || z <= Z_buf(x, y))
        return false;

    image(x, y) = PaintColor;
	Z_buf(x, y) = z; 

    return true;
}

bool Graphics::drawPoint(Mat<ARGB>& image, Mat<int>& Z_buf, int x, int y, int z, double nx, double ny, double nz)
{
    if (image.isOut(x, y) || z <= Z_buf(x, y))
        return false;

    double a = Illumination::Phong(nx, ny, nz);

    image(x, y) = RGB::mul(PaintColor, a);;
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
		delta[d] *= inc[d];
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

void Graphics::drawTriangle(Mat<ARGB>& image, Mat<int>& Z_buf, vector<int>& p1, vector<int>& p2, vector<int>& p3) {
    vector<int>* q[3] = { &p1, &p2, &p3 };

    if ((*q[0])[1] > (*q[1])[1]) swap(q[0], q[1]);
    if ((*q[0])[1] > (*q[2])[1]) swap(q[0], q[2]);
    if ((*q[1])[1] > (*q[2])[1]) swap(q[1], q[2]);

    vector<int> delta[3], a[3], inc[3], p[3];
    {
        p[0] = *(q[0]);
        p[1] = *(q[0]);
        p[2] = *(q[1]);
        
        for (int i = 0; i < 3; i++) {
            delta[i].resize(3);
            a[i].resize(3);
            inc[i].resize(3);
            fill(a[i].begin(), a[i].end(), 0);
        }

        for (int d = 0; d < 3; d++) {
            delta[0][d] = (*q[1])[d] - (*q[0])[d];
            delta[1][d] = (*q[2])[d] - (*q[0])[d];
            delta[2][d] = (*q[2])[d] - (*q[1])[d];

            for (int i = 0; i < 3; i++) {
                inc[i][d] = delta[i][d] == 0 ? 0 : (delta[i][d] > 0 ? 1 : -1);
                delta[i][d] *= inc[i][d];
            }
        }
    }

    ARGB color = PaintColor;
    {
        vector<double> d1(3), d2(3);
        for (int d = 0; d < 3; d++) {
            d1[d] = p3[d] - p1[d];
            d2[d] = p2[d] - p1[d];
        }
        double a = Illumination::Phong(
            d1[1] * d2[2] - d1[2] * d2[1],
            d1[2] * d2[0] - d1[0] * d2[2],
            d1[0] * d2[1] - d1[1] * d2[0]);

        PaintColor = RGB::mul(PaintColor, a);
    }

    for (int l = 0; l <= 2; l += 2) {
        if (delta[l][1] == 0)
            continue;

        for (int i = (*q[l == 0 ? 0 : 1])[1]; i < (*q[l == 0 ? 1 : 2])[1]; i++) {
            for (int d = 0; d < 3; d++) {
                a[l][d] += delta[l][d];
                while (a[l][d] >= delta[l][1]) {
                    a[l][d] -= delta[l][1];
                    p[l][d] += inc[l][d];
                }

                a[1][d] += delta[1][d];
                while (a[1][d] >= delta[1][1]) {
                    a[1][d] -= delta[1][1];
                    p[1][d] += inc[1][d];
                }
            }
            if (p[1][0] != p[l][0]) {
                int st = p[l][0] >= p[1][0] ? 1 : l,
                    ed = p[l][0] >= p[1][0] ? l : 1,
                    dx = abs(p[l][0] - p[1][0]),
                    dz = (p[l][2] >= p[1][2]) ? 1 : -1,
                    k = p[1][2], e = 0;

                for (int j = p[st][0]; j <= p[ed][0]; j++) {
                    e += dz * (p[l][2] - p[1][2]);

                    while (e >= dx) {
                        e -= dx;
                        k += dz;
                    }
                    drawPoint(image, Z_buf, j, i, k);
                }
            }
        }
    }
    PaintColor = color;
}


void Graphics::drawTriangle(Mat<ARGB>& image, Mat<int>& Z_buf, vector<double>& p1, vector<double>& p2, vector<double>& p3) {
    vector<int> pi1(3), pi2(3), pi3(3);

    pi1 = { (int)p1[0], (int)p1[1], (int)p1[2] };
    pi2 = { (int)p2[0], (int)p2[1], (int)p2[2] };
    pi3 = { (int)p3[0], (int)p3[1], (int)p3[2] };

    vector<int>* q[3] = { &pi1, &pi2, &pi3 };

    if ((*q[0])[1] > (*q[1])[1]) swap(q[0], q[1]);
    if ((*q[0])[1] > (*q[2])[1]) swap(q[0], q[2]);
    if ((*q[1])[1] > (*q[2])[1]) swap(q[1], q[2]);

    vector<int> delta[3], a[3], inc[3], p[3];
    {
        p[0] = *(q[0]);
        p[1] = *(q[0]);
        p[2] = *(q[1]);

        for (int i = 0; i < 3; i++) {
            delta[i].resize(3);
            a[i].resize(3);
            inc[i].resize(3);
            fill(a[i].begin(), a[i].end(), 0);
        }

        for (int d = 0; d < 3; d++) {
            delta[0][d] = (*q[1])[d] - (*q[0])[d];
            delta[1][d] = (*q[2])[d] - (*q[0])[d];
            delta[2][d] = (*q[2])[d] - (*q[1])[d];

            for (int i = 0; i < 3; i++) {
                inc[i][d] = delta[i][d] == 0 ? 0 : (delta[i][d] > 0 ? 1 : -1);
                delta[i][d] *= inc[i][d];
            }
        }
    }

    ARGB color = PaintColor;
    {
        vector<double> d1(3), d2(3);
        Matrix::sub(d1, p3, p1);
        Matrix::sub(d2, p2, p1);

        double a = Illumination::Phong(
            d1[1] * d2[2] - d1[2] * d2[1],
            d1[2] * d2[0] - d1[0] * d2[2],
            d1[0] * d2[1] - d1[1] * d2[0]);

        PaintColor = RGB::mul(PaintColor, a);
    }

    for (int l = 0; l <= 2; l += 2) {
        if (delta[l][1] == 0)
            continue;

        for (int i = (*q[l == 0 ? 0 : 1])[1]; i < (*q[l == 0 ? 1 : 2])[1]; i++) {
            for (int d = 0; d < 3; d++) {
                a[l][d] += delta[l][d];
                while (a[l][d] >= delta[l][1]) {
                    a[l][d] -= delta[l][1];
                    p[l][d] += inc[l][d];
                }

                a[1][d] += delta[1][d];
                while (a[1][d] >= delta[1][1]) {
                    a[1][d] -= delta[1][1];
                    p[1][d] += inc[1][d];
                }
            }
            if (p[1][0] != p[l][0]) {
                int st = p[l][0] >= p[1][0] ? 1 : l,
                    ed = p[l][0] >= p[1][0] ? l : 1,
                    dx = abs(p[l][0] - p[1][0]),
                    dz = (p[l][2] >= p[1][2]) ? 1 : -1,
                    k = p[1][2], e = 0;

                for (int j = p[st][0]; j <= p[ed][0]; j++) {
                    e += dz * (p[l][2] - p[1][2]);

                    while (e >= dx) {
                        e -= dx;
                        k += dz;
                    }
                    drawPoint(image, Z_buf, j, i, k);
                }
            }
        }
    }
    PaintColor = color;
}

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