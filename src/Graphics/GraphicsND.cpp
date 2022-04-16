#include "GraphicsND.h"

// ����
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
*                    �������� any-D
*	[�㷨]: ���ö����Ʊ�ʾ�����������λ�ã�
			��С����Ϊȫ0����󶥵�Ϊȫ1
			0: ��ά������ֵ == ��С�����Ӧֵ
			1: ��ά������ֵ == ��󶥵��Ӧֵ
*	[����]:
		[1] �Զ�����˳��������ж���
			[2] ���Ӹõ�����бȸõ�����1�ĵ�
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

/*--------------------------------[ ������ ]--------------------------------
*	[����]:
		[1] ����ÿһ����������
		[2] ���Ƹø���Ӧ��, ��ά�ȷ���Ĵ�min[dim] -> max[dim]��ֱ�߶�
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