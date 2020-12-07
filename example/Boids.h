#ifndef BOIDS_H
#define BOIDS_H
#include "LiGu_Graphics/Graphics3D.h"
#include "LiGu_AlgorithmLib/ComputationalGeometry.h"
#include <ctime> 
void Delay(int time){clock_t now = clock(); while (clock() - now < time); }//time*1000为秒数 
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
	Graphics3D g{ 1000, 1000 };
	double size = 300;
	double visualField = 20, visualAngle = -0.5;		// 能见范围 // 能见角度cosθ
	double deltaT = 1, speed = 3;
	double weight[4] = { 1,1,1,4 };					// 各规则的权值
	Mat<double>* obstacleAvoidance_direction;
	/*********************************[ 函数 ]*********************************/
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
		int numViewDirections;  Mat<double>zero(3,1);
		obstacleAvoidance_direction = getSphereFibonacciPoint(numViewDirections);
	}
	void obstacleAvoidance(Boids_cell& cell) {
		// 确定坐标变换矩阵
		Mat<double>rorateAxis(3, 1), rotateMat, Temp(3, 1);
		rorateAxis.crossProduct(cell.v, obstacleAvoidance_direction[0]);
		double theta = -acos(Temp.dot(cell.v, obstacleAvoidance_direction[0])); // cosθ = a * b / (|a| * |b|)
		g.rotate(rorateAxis.normalization(), theta, Temp, rotateMat);

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
		for (int i = 0; i < 3; i++) {
			if (fabs(point[i]) >= size)return true;
		}
		return false;
	}
	/*--------------------------------[ Boidss 主函数 ]--------------------------------*/
	void Boids_main() {
		// set Graphics3D
		Mat<double> Axis(3, 1); Axis[0] = 0.7; Axis[1] = 0.4; Axis[2] = 0.59;
		Mat<double> center(3, 1); center[0] = 500; center[1] = 500; center[2] = 500;
		double theta = 1;
		// 
		Mat<double> zero(3, 1), CuboidL(3, 1), CuboidR(3, 1);
		CuboidR[0] = size; CuboidR[1] = size; CuboidR[2] = size; CuboidL[0] = -size; CuboidL[1] = -size; CuboidL[2] = -size;
		// set Boids_cell
		int N = 200;
		Boids_cell cell[200];
		for (int i = 0; i < N; i++) {
			cell[i].r.rands(3, 1, -size / 2, size / 2);
			cell[i].v.rands(3, 1, 0, 1);
			cell[i].v.normalization();
		}
		Mat<double> temp;
		obstacleAvoidance_init();
		while (1) {
			theta += 0.01; if (theta >= 2 * 3.1415926) theta = 0;
			g.TransformMat.E(4);
			g.rotate(Axis, theta, zero);
			g.translation(center);

			play(cell, N);
			for (int i = 0; i < N; i++)
				g.drawLine(cell[i].r, temp.add(temp.mult(10, cell[i].v), cell[i].r));

			Delay(20);
			g.g->PicWrite("D:/LIGU.ppm");
			g.g->clear(0);
			g.drawCuboid(CuboidL, CuboidR);
		}
	}
};



/*
int main()
{
	Delay(1000);
	Boids boids;
	boids.Boids_main();
}
*/
#endif
