#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "LiGu_AlgorithmLib/Mat.h"
#define PI 3.141592653589
FILE* fin;
/*-------------------------------- Perlin Noise --------------------------------
Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
--------------------------------------------------------------------------------*/
double PerlinNoise(double x, double y, Mat<double>& randomGridGradient) {
	// 对四个格点
	int x0[] = { x, x + 1, x, x + 1 }, y0[] = { y, y, y + 1, y + 1 };
	double n[4];
	for (int i = 0; i < 4; i++) {
		//[1] 格点随机梯度矢量
		double random = randomGridGradient(x0[i], y0[i]);
		//[2] (x,y)与格点距离,梯度点积
		double dx = x - x0[i], dy = y - y0[i];
		n[i] = dx * cos(random) + dy * sin(random);
	}
	//[3] 插值
	double sx = x - (int)x, sy = y - (int)y;
	double ix0 = (n[1] - n[0]) * (3.0 - sx * 2.0) * sx * sx + n[0];
	double ix1 = (n[3] - n[2]) * (3.0 - sx * 2.0) * sx * sx + n[2];
	return (ix1 - ix0) * (3.0 - sy * 2.0) * sy * sy + ix0;
}
void PerlinNoise(Mat<double>& map, int frequency) {
	Mat<double> randomGridGradient;
	randomGridGradient.rands(frequency + 1, frequency + 1, 0, 256);
	for (int y = 0; y < map.rows; y++)
		for (int x = 0; x < map.cols; x++)
			map(x, y) = PerlinNoise((double)x / map.cols * frequency, (double)y / map.rows * frequency, randomGridGradient);
}
/*-------------------------------- Mountain --------------------------------
Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
--------------------------------------------------------------------------------*/
void Mountain(Mat<double>& map) {
	Mat<double> t(map.rows, map.cols);
	for (int k = 2; k <= map.rows; k *= 2) {
		PerlinNoise(t, k);
		map.add(map, t.mult(1.0 / k, t));
	}map.mult(map.rows, map);
}
/*-------------------------------- translation --------------------------------*/
void translation(Mat<double>& delta, Mat<double>& transMat) {
	int n = delta.rows;
	Mat<double> translationMat(n + 1);
	for (int i = 0; i < n; i++)translationMat(i, n) = delta[i];
	transMat.mult(translationMat, transMat);
}
/*-------------------------------- rotate --------------------------------
*	[公式]: v' = q v q`¹
		q = [cos(θ/2), u s in(θ/2)]
		v=>[0,v]经旋转轴u旋转Ѳ角后得到v'
	多次旋转:
		v' = q1q2 v q2`¹q1`¹ = (q1q2) v (q1q2)`¹
	四元数化旋转矩阵:
		四元数左乘:
		q v =	[a -b -c -d] v
				|b  a -d  c|
				|c  d  a -b|
				[d -c  b  a]
		四元数右乘:
		v q =	[a -b -c -d] v
				|b  a  d -c|
				|c -d  a  b|
				[d  c -b  a]
四元数乘法: q1q2 = (v1 × v2 + w1v2 + w2v1 , w1w2 - v1·v2)
----------------------------------------------------------------*/
void rotate(Mat<double>& rotateAxis, double rotateAngle, Mat<double>& center, Mat<double>& transMat) {
	Mat<double> tmp(3, 1);
	if (rotateAxis == tmp) return;
	translation(center.negative(tmp), transMat);
	// rotate
	rotateAxis.normalization();
	Mat<double> q(4, 1);
	{ double t[] = { cos(rotateAngle / 2)
		,rotateAxis[0] * sin(rotateAngle / 2)
		,rotateAxis[1] * sin(rotateAngle / 2)
		,rotateAxis[2] * sin(rotateAngle / 2) }; q.getData(t); }
	// rotate Mat
	Mat<double> rotateMat(4);
	{double t[] = {
		1 - 2 * q[2] * q[2] - 2 * q[3] * q[3],2 * q[1] * q[2] - 2 * q[0] * q[3],2 * q[1] * q[3] + 2 * q[0] * q[2],0,
		2 * q[1] * q[2] + 2 * q[0] * q[3],1 - 2 * q[1] * q[1] - 2 * q[3] * q[3],2 * q[2] * q[3] - 2 * q[0] * q[1],0,
		2 * q[1] * q[3] - 2 * q[0] * q[2], 2 * q[2] * q[3] + 2 * q[0] * q[1],1 - 2 * q[1] * q[1] - 2 * q[2] * q[2],0
	}; rotateMat.getData(t); }
	transMat.mult(rotateMat, transMat);
	translation(center, transMat);
}
/*-------------------------------- Fractal Tree 3D --------------------------------*/
void FractalTree3D(std::vector<Mat<double>>& linesSt, std::vector<Mat<double>>& linesEd, int level, double alpha) {
	if (level <= 0)return;
	Mat<double> lineSt0 = linesSt.back();
	Mat<double> lineEd0 = linesEd.back();
	// 确定旋转矩阵
	Mat<double> direction, rotateAxis, rotateMat(4), Temp(3, 1);
	direction.add(lineEd0, lineSt0.negative(direction));
	{ double t[] = { 0, 0, 1 }; Temp.getData(t); }
	rotateAxis.crossProduct(direction, Temp);
	double rotateAngle = -acos(Temp.dot(direction, Temp) / direction.norm()); // cosθ = a * b / (|a| * |b|)
	rotate(rotateAxis, rotateAngle, Temp.zero(3, 1), rotateMat);
	Temp = rotateMat; rotateMat.E(3);
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)rotateMat(i, j) = Temp(i, j);
	//递归
	int fork = 3;
	double Lenth = direction.norm();
	Mat<double> endPoint(3, 1);
	for (int i = 0; i < fork; i++) {
		{double t[] = { sin(alpha) * cos((double)i * 2 * PI / fork), sin(alpha) * sin((double)i * 2 * PI / fork), cos(alpha) }; endPoint.getData(t); }
		endPoint.add(lineEd0, endPoint.mult(0.7 * Lenth, endPoint.mult(rotateMat, endPoint)));
		linesSt.push_back(lineEd0);
		linesEd.push_back(endPoint);
		FractalTree3D(linesSt, linesEd, level - 1, alpha);
	}
}
/*-------------------------------- obj Save --------------------------------
Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
--------------------------------------------------------------------------------*/
void ObjSave(Mat<double>& map, const char* url) {
	static int index = 1;
	for (int y = 0; y < map.rows; y++)
		for (int x = 0; x < map.cols; x++)
			fprintf(fin, "v %f %f %f\n", (float)x, map(x, y), (float)y);
	for (int y = 1; y < map.rows; y++)
		for (int x = 1; x < map.cols; x++)
			fprintf(fin, "f %d %d %d %d\n", (y - 1) * map.cols + x + index, (y - 1) * map.cols + x - 1 + index, y * map.cols + x - 1 + index, y * map.cols + x + index);
	index += map.rows * map.cols;
}



/*-------------------------------- draw Triangle --------------------------------
Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
--------------------------------------------------------------------------------*/
void drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3) {
	fprintf(fin, "v %f %f %f\n", p1[0], p1[2], p1[1]);
	fprintf(fin, "v %f %f %f\n", p2[0], p2[2], p2[1]);
	fprintf(fin, "v %f %f %f\n", p3[0], p3[2], p3[1]);
	fprintf(fin, "f -3 -2 -1\n");
}
/*-------------------------------- draw Rectangle --------------------------------
Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
--------------------------------------------------------------------------------*/
void drawRectangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4) {
	drawTriangle(p1, p2, p3); drawTriangle(p1, p3, p4);
}
/*-------------------------------- draw cylinder --------------------------------
Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
--------------------------------------------------------------------------------*/
void drawCylinder(Mat<double>& st, Mat<double>& ed, double r) {
	Mat<double> direction, rotateAxis(3, 1), rotateMat(4), tmp(3, 1);
	direction.add(st, ed.negative(direction));
	{ double t[] = { 0, 0, 1 }; tmp.getData(t); }
	rotateAxis.crossProduct(direction, tmp);
	double rotateAngle = -acos(tmp.dot(direction, tmp) / direction.norm());
	rotate(rotateAxis, rotateAngle, tmp.zero(3, 1), rotateMat);
	tmp = rotateMat; rotateMat.E(3);
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)rotateMat(i, j) = tmp(i, j);
	//
	Mat<double> stPoint(3, 1), preStPoint;
	Mat<double> edPoint(3, 1), preEdPoint;
	for (int i = 0; i <= 360; i++) {
		Mat<double> delta(3, 1);
		{double t[] = { r * cos((double)i * 2 * PI / 360), r * sin((double)i * 2 * PI / 360),0 }; delta.getData(t); }
		delta.mult(rotateMat, delta);
		stPoint.add(st, delta);
		edPoint.add(ed, delta);
		if (i != 0) {
			drawTriangle(st, stPoint, preStPoint);
			drawTriangle(ed, edPoint, preEdPoint);
			drawRectangle(stPoint, edPoint, preEdPoint, preStPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}
void drawCylinder(Mat<double>& st, Mat<double>& ed, double stR, double edR) {
	Mat<double> direction, rotateAxis(3, 1), rotateMat(4), tmp(3, 1);
	direction.add(st, ed.negative(direction));
	{double t[] = { 0, 0, 1 }; tmp.getData(t); }
	rotateAxis.crossProduct(direction, tmp);
	double rotateAngle = -acos(tmp.dot(direction, tmp));
	rotate(rotateAxis, rotateAngle, tmp.zero(3, 1), rotateMat);
	tmp = rotateMat; rotateMat.E(3);
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)rotateMat(i, j) = tmp(i, j);
	//
	Mat<double> stPoint(3, 1), preStPoint;
	Mat<double> edPoint(3, 1), preEdPoint;
	for (int i = 0; i <= 360; i++) {
		Mat<double> delta(3, 1);
		{double t[] = { cos((double)i * 2 * PI / 360), sin((double)i * 2 * PI / 360),0 }; delta.getData(t); }
		delta.mult(rotateMat, delta);
		stPoint.add(st, stPoint.mult(stR, delta));
		edPoint.add(ed, edPoint.mult(edR, delta));
		if (i != 0) {
			drawTriangle(st, stPoint, preStPoint);
			drawTriangle(ed, edPoint, preEdPoint);
			drawRectangle(stPoint, edPoint, preEdPoint, preStPoint);
		}
		preStPoint = stPoint;
		preEdPoint = edPoint;
	}
}


#include<time.h>
int main() {
	time_t timenow;
	srand((unsigned)time(&timenow));
	fin = fopen("D:/LIGU.obj", "w+");
	// Mountain
	Mat<double> mountain(1024, 1024);
	Mountain(mountain);
	ObjSave(mountain, "D:/LIGU.obj");
	// Sea
	Mat<double> sea(1024, 1024);
	ObjSave(sea, "D:/LIGU.obj");
	// Tree
	std::vector<Mat<double>> TreeSt, TreeEd;
	{
		int indexMax;
		mountain.max(indexMax);
		Mat<double> tMat(3, 1);
		{double t[] = { indexMax / mountain.rows,indexMax % mountain.rows,mountain[indexMax] }; tMat.getData(t); }
		TreeSt.push_back(tMat);
		{double t[] = { indexMax / mountain.rows,indexMax % mountain.rows,mountain[indexMax] + 1.5 }; tMat.getData(t); }
		TreeEd.push_back(tMat);
	}
	FractalTree3D(TreeSt, TreeEd, 6, (double)30 * 2 * PI / 360);
	for (int i = 0; i < TreeSt.size(); i++) {
		Mat<double> tmp;
		double Width = 0.1 * (tmp.add(TreeSt[i], TreeEd[i].negative(tmp))).norm();
		drawCylinder(TreeSt[i], TreeEd[i], Width, 0.7 * Width);
	}
}




/******************************************************************************
*                    Boids
*	[Rule]:
		[1] Collision Avoidance: avoid collisions with nearby flockmates
		[2] Velocity Matching: attempt to match velocity with nearby flockmates
		[3] Flock Centering: attempt to stay close to nearby flockmates
*	[过程]:
		[可见域内]:
			[1] 可见距离判定
			[2] 可见角度判定
		[Rule 1]: 可见域内，同各个个体的反向方向，除以距离，作为加速度 a1
		[Rule 2]: 可见域内, 质心速度, 作为加速度 a2
		[Rule 3]: 可见域内, 质心位置方向, 作为加速度 a3
		[障碍规避]:
			[1] 生成均匀球面点分布
			[2] 依次从速度方向,向四周,发出射线进行障碍检测,射线距离在可见域内
			[3] 选择第一束探测无障碍的射线, 作为避障加速度 a4
		[ v ]: v = (Σwi·ai)·dt
				v 归一化单位矢量: ^v  =>  v = v0·^v
		[ r ]: r = v·dt = v0·dt·^v
*	[Referance]:
		[1] Thanks for https://github.com/SebLague/Boids
******************************************************************************/
class Boids {
public:
	struct Boids_cell { Mat<double> r{ 3,1 }, v{ 3,1 }, a{ 3,1 }; };
	double size = 300, deltaT = 1, speed = 3;
	double visualField = 20, visualAngle = -0.5;		// 能见范围 // 能见角度cosθ
	double weight[4] = { 1,1,1,4 };					// 各规则的权值
	Mat<double>* obstacleAvoidance_direction;
	/*********************************[ 函数 ]*********************************/
	/*--------------------------------[ Boids ]--------------------------------*/
	Boids() {
		int N = 200;
		Boids_cell cell[200];
		for (int i = 0; i < N; i++) {
			cell[i].r.rands(3, 1, -size / 2, size / 2);
			cell[i].v.rands(3, 1, 0, 1);
			cell[i].v.normalization();
		}
		obstacleAvoidance_init();
	}
	/*--------------------------------[ play 运行 ]--------------------------------*/
	void play(Boids_cell cell[], int n) {
		for (int i = 0; i < n; i++) {
			rule(cell, n, i);
			obstacleAvoidance(cell[i]);
		}
		Mat<double> temp(3, 1);
		for (int i = 0; i < n; i++) {
			cell[i].v.add(cell[i].v, temp.mult(deltaT, cell[i].a));
			cell[i].v.normalization();
			cell[i].r.add(cell[i].r, temp.mult(deltaT * speed, cell[i].v));
		}
	}
	/*--------------------------------[ rule 规则 ]--------------------------------*/
	void rule(Boids_cell cell[], int n, int index) {
		Mat<double> distance(3, 1);
		Mat<double> avoidDirection(3, 1), groupVelocity(3, 1), groupCenter(3, 1), temp;
		int groupNum = 0;
		for (int i = 0; i < n; i++) {
			if (i == index)continue;
			// 能见范围
			distance.add(cell[i].r, cell[index].r.negative(distance));
			if (distance.norm() > visualField)continue;
			// 能见角度
			double angle = distance.dot(distance, cell[index].v) / (distance.norm() * cell[index].v.norm()); // cosθ = a * b / (|a| * |b|)
			if (angle < visualAngle)continue;
			groupNum++;
			//[Rule 1]: collisionAvoid  [Rule 2]: velocityMatching  [Rule 3]: flockCentering
			avoidDirection.add(avoidDirection, temp.mult(1 / distance.norm(), distance));
			groupVelocity.add(groupVelocity, cell[i].v);
			groupCenter.add(groupCenter, distance);
		}
		avoidDirection.negative(avoidDirection);
		// Update Acceleration
		cell[index].a.clean();
		if (groupNum == 0) return;
		cell[index].a.add(cell[index].a, avoidDirection.mult(weight[0], avoidDirection.normalization()));
		cell[index].a.add(cell[index].a, groupVelocity.mult(weight[1], groupVelocity.normalization()));
		cell[index].a.add(cell[index].a, groupCenter.mult(weight[2], groupCenter.normalization()));
	}
	/*--------------------------------[ obstacleAvoidance 障碍规避 ]--------------------------------*/
	void obstacleAvoidance_init() {// 均匀球面点
		int numViewDirections;  Mat<double>zero(3, 1);
		//obstacleAvoidance_direction = getSphereFibonacciPoint(numViewDirections);
	}
	void obstacleAvoidance(Boids_cell& cell) {
		// 确定坐标变换矩阵
		Mat<double>rorateAxis(3, 1), rotateMat, Temp(3, 1);
		rorateAxis.crossProduct(cell.v, obstacleAvoidance_direction[0]);
		double theta = -acos(Temp.dot(cell.v, obstacleAvoidance_direction[0])); // cosθ = a * b / (|a| * |b|)
		rotate(rorateAxis.normalization(), theta, Temp, rotateMat);

		Temp = rotateMat; rotateMat.E(3);
		for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)rotateMat(i, j) = Temp(i, j);
		// 碰撞检测
		Mat<double>derection(3, 1);
		for (int k = 0; k < 200; k++) {
			derection.mult(rotateMat, obstacleAvoidance_direction[k]);
			bool flag = 1; int distance = visualField;
			while (--distance > 0) {
				Temp.add(Temp.mult(distance, derection), cell.r);
				if (obstacleFunction(Temp)) { flag = 0; break; }
			}
			if (flag)break;
		}
		cell.a.add(cell.a, derection.mult(weight[3], derection));
	}
	bool obstacleFunction(Mat<double>& point) {
		for (int i = 0; i < 3; i++)
			if (fabs(point[i]) >= size)return true;
		return false;
	}
};