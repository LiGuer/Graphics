#include "GraphicsND.h"

// 画线
void Graphics::drawPoint(Mat<ARGB>& image, Mat<>& p0) {
	static Mat<> point;

	point.zero(TransformMat.rows);
	point[0] = 1;
	for (int i = 0; i < p0.rows; i++)
		point[i + 1] = p0[i];

	Matrix::mul(point, TransformMat, point);

	drawPoint(image, (int)point(1), (int)point(2));
}

void Graphics::drawPoint(Mat<ARGB>& image, double x, double y, double z) {
	static Mat<> point;

	point.zero(TransformMat.rows);
	point[0] = 1;
	if (x != 0) point[1] = x;
	if (y != 0) point[2] = y;
	if (z != 0) point[3] = z;

	Matrix::mul(point, TransformMat, point);

	drawPoint(image, (int)point(1), (int)point(2));
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
void Graphics::drawSuperCuboid(Mat<ARGB>& image, Mat<>& pMin, Mat<>& pMax) {
	unsigned int Dim = pMin.rows, maxCode = 0;
	Mat<> st, ed;
	for (int i = 0; i < Dim; i++) 
		maxCode += 1 << i;

	for (unsigned int code = 0; code < maxCode; code++) {
		st = pMin;
		for (int i = 0; i < Dim; i++)
			if (code & (1 << i))
				st[i] = pMax[i];

		ed = st;
		for (int i = 0; i < Dim; i++) {
			if (ed[i] == pMin[i]) {
				ed[i] = pMax[i];
				//drawLine(image, st, ed);
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
void Graphics::drawGrid(Mat<ARGB>& image, Mat<>& delta, Mat<>& max, Mat<>& min) {
	int times = 1, cur;
	for (int dim = 0; dim < min.rows; dim++) 
		times *= (max[dim] - min[dim]) / delta[dim] + 1;

	Mat<> point(min), st, ed; point[0] -= delta[0];
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
			for (int dim = 0; dim < min.rows; dim++) {
				st[dim] = min[dim];
				ed[dim] = max[dim];
				//drawLine(image, st, ed);
				st[dim] = point[dim];
				ed[dim] = point[dim];
			}
		}
	}
}