#include "GraphicsND.h"

/*
 * 画点
 */ 
bool Graphics::drawPoint (Mat<ARGB>& image, Mat<int>* Z_buf, vector<int>& p)
{
	int dim = p.size();

    if (image.isOut(p[0], p[1]))
        return false;

	for (int d = 2; d < dim; d++)
	    if (p[d] < Z_buf[d-2](p[0], p[1]))
        	return false;

    image(p[0], p[1]) = PaintColor;

	for (int d = 2; d < dim; d++)
	    Z_buf[d-2] = p[d]; 

    return true;
}

void Graphics::drawPoint(Mat<ARGB>& image, vector<double>& p0) {
	static vector<double> point;

	point.resize(TransformMat.rows, 0);
	point[0] = 1;
	for (int i = 0; i < p0.size(); i++)
		point[i + 1] = p0[i];

	Matrix::mul(point, TransformMat, point);

	drawPoint(image,
		(int)(image.rows / 2 + point[1]),
		(int)(image.cols / 2 + point[2])
	);
}

/*
 * 画线
 */

void Graphics::drawLine (Mat<ARGB>& image, Mat<int>* Z_buf, vector<double>& st, vector<double>& ed)
{
	int dim = st.size();
	vector<double>
		a(dim),
		inc(dim),
		delta(dim),
		p = st;
	vector<int> p_int(dim);

	Matrix::sub(delta, ed, st);
 
	for (int d = 0; d < dim; d++) {
		inc[d] = delta[d] == 0 ? 0 : (delta[d] > 0 ? 1 : -1);
		delta[d] *= inc[d];    // |Δ| 
	}
	int max_delta = Matrix::max(delta);

	for (int i = 0; i < p.size(); i++)
		p_int[i] = (int) p[i];

	for (int i = 0; i <= max_delta; i++) {
		drawPoint(image, Z_buf, p_int);

		for (int d = 0; d < dim; d++) {
			a[d] += delta[d];
			if (a[d] >= max_delta) {
				a[d] -= max_delta;
				p[d] += inc[d];
			}
		}
	}
}


void Graphics::drawLine(Mat<ARGB>& image, vector<double>& st, vector<double>& ed) {
	static vector<double> p1, p2;

	p1.resize(TransformMat.rows, 0);
	p2.resize(TransformMat.rows, 0);

	p1[0] = 1;
	for (int i = 0; i < st.size(); i++)
		p1[i + 1] = st[i];

	p2[0] = 1;
	for (int i = 0; i < ed.size(); i++)
		p2[i + 1] = ed[i];

	Matrix::mul(p1, TransformMat, p1);
	Matrix::mul(p2, TransformMat, p2);

	drawLine(image, 
		image.rows / 2 + p1[1], image.rows / 2 + p1[2],
		image.cols / 2 + p2[1], image.cols / 2 + p2[2]
	);
}

/*----------------------------------------------------------------
*                    画立方体 any-D
*	[算法]: 利用二进制表示各顶点的坐标位置，
			最小顶点为全0，最大顶点为全1
			0: 该维度坐标值 == 最小顶点对应值
			1: 该维度坐标值 == 最大顶点对应值
*	[流程]:
		[1] 以二进制顺序遍历所有顶点
			[2] 连接该点和所有比该点编码多1的点
----------------------------------------------------------------*/
void Graphics::drawSuperCuboid(Mat<ARGB>& image, vector<double>& pMin, vector<double>& pMax) {
	unsigned int Dim = pMin.size(), maxCode = (1 << Dim) - 1;
	vector<double> st, ed;

	for (unsigned int code = 0; code < maxCode; code++) {
		st = pMin;
		for (int i = 0; i < Dim; i++)
			if (code & (1 << i))
				st[i] = pMax[i];

		ed = st;
		for (int i = 0; i < Dim; i++) {
			if (ed[i] == pMin[i]) {
				ed[i] = pMax[i];
				drawLine(image, st, ed);
				ed[i] = pMin[i];
			}
		}
	}
}


/*--------------------------------[ 画网格 ]--------------------------------
*	[过程]:
		[1] 计算每一个格点的坐标
		[2] 绘制该格点对应的, 各维度方向的从min[dim] -> max[dim]的直线段
---------------------------------------------------------------------------*/
void Graphics::drawGrid(Mat<ARGB>& image, vector<double>& delta, vector<double>& max, vector<double>& min) {
	int times = 1, cur;
	for (int dim = 0; dim < min.size(); dim++) 
		times *= (max[dim] - min[dim]) / delta[dim] + 1;

	vector<double> point(min), st, ed; point[0] -= delta[0];
	for (int i = 0; i < times; i++) {
		//[1]
		cur = 0; point[cur] += delta[cur];
		while (point[cur] > max[cur]) {
			point[cur] = min[cur]; cur++;
			point[cur] += delta[cur];
		}

		//[2]
		//if (!LINE)
			//drawPoint(image, point);

		//if (LINE)
		{
			st = point;
			ed = point;
			for (int dim = 0; dim < min.size(); dim++) {
				st[dim] = min[dim];
				ed[dim] = max[dim];
				drawLine(image, st, ed);
				st[dim] = point[dim];
				ed[dim] = point[dim];
			}
		}
	}
}