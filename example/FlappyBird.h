#include <ctime> 
void Delay(int time) { clock_t now = clock(); while (clock() - now < time); }//time*1000为秒数 

#include "LiGu_Graphics/Graphics.h"

/*--------------------------------[ QLearning ]--------------------------------
*	[定义]:Q(s,a) = (1 + lr)·Q(s,a) + lr·( R + g·max Q(s',:) )
			s: state	a: action	R: reward	lr: learning rate	g: forget factor
*	[原理]:
		选择动作: ε-greedy方法: 
			每个状态以ε概率随机选取动作，1-ε概率选择当前最优解
		眼前利益R: 反馈值.
		记忆中的利益 max Q(s',:): 小鸟记忆里，新位置[公式]能给出的最大效用值.
		forget factor越大，小鸟就会越重视以往经验，越小，小鸟只重视眼前利益R.
*	[流程]:
		[1] Init Q table arbitrarily
		[2] Repeat (for each episode), until s is terminal
			[3] Choose a from s using policy derived from Q (eg. ε-greedy)
			[4] Take action a, observe r s'
			[5] Q(s,a) = (1 + lr)·Q(s,a) + lr·( R + g·max Q(s',:) )
				s = s'
*	[Ps]:
		可以逐渐降低随机选取动作的概率ε，一开始随机率可达100%
			然后随训练次数的深入，应当逐渐降低随机概率。
-----------------------------------------------------------------------------*/
class QLearning {
public:
	double learnRate = 0.6, Gamma = 0.8, greedy = 0.9; //奖励递减值# 贪婪度
	int actionNum = 0, stateNum = 0;
	Mat<double> QTable;
	double preState = 0;
	/*---------------- 初始化 ----------------*/
	QLearning(int _stateNum, int _actionNum) { init(_stateNum, _actionNum); }
	void init(int _stateNum, int _actionNum) {
		actionNum = _actionNum;
		stateNum = _stateNum;
		QTable.zero(_stateNum, _actionNum);
	}
	/*---------------- 选择行为 ----------------*/
	int chooseAction(int state) {
		int action = 0;
		bool flag = 1;
		for (int i = 0; i < actionNum; i++)
			if (QTable(state, i) != 0) { flag = 0; break; }
		if (rand() / double(RAND_MAX) < greedy || flag) return rand() % actionNum;
		double t = -DBL_MAX;
		for (int i = 0; i < actionNum; i++)
			if (QTable(state, i) > t) { t = QTable(state, i); action = i; }
		return action;
	}
	/*---------------- 反馈学习 ----------------*/
	void feedback(int state, int action, double R) {
		double t = -DBL_MAX;
		for (int i = 0; i < actionNum; i++)
			t = QTable(state, i) > t ? QTable(state, i) : t;
		QTable(preState, action) = (1 - learnRate) * QTable(preState, action) + learnRate * (R + Gamma * t);
		preState = state;
	}
};
/*--------------------------------[ FlappyBird ]--------------------------------
*	[状态空间]: [-xb, x0-xb], [-y0, y0]
-----------------------------------------------------------------------------*/
class FlappyBird {
public:
	int windowWidth = 400, windowHeight = 600;
	Graphics g{ windowWidth, windowHeight };
	int birdY = windowHeight / 2, birdX = 100, pillarsWidth = 80, holeHeight = 100, pillarsX = windowWidth, pillarsY = 0;
	int pillarsSpeed = 10, BirdSpeed = 20;
	QLearning q{ (windowWidth + pillarsWidth) / (2 * pillarsSpeed) * (windowHeight * 2) / BirdSpeed, 2 };
	void init() {
		birdY = windowHeight / 2; initPillars();
	}
	void initPillars() {
		pillarsX = windowWidth;
		pillarsY = rand() / double(RAND_MAX) * (windowHeight - 100 - 100) + 100;	//[st,ed)
	}
	void play() {
		init();
		bool flag = 0;
		long long time = 0;
		while (true) {
			//pillars &  boid
			pillarsX -= pillarsSpeed;
			birdY -= BirdSpeed;
			if (flag) {
				int dx = pillarsX - birdX + birdX + pillarsWidth;
				int dy = pillarsY - birdY + windowHeight;
				int state = dy / BirdSpeed * windowWidth / (2 * pillarsSpeed) + dx / (2 * pillarsSpeed);
				//printf("%d %d %d %d %d\n", pillarsY, birdY,dx, dy, state);
				int action = q.chooseAction(state);
				if (action)birdY += 3 * BirdSpeed;
				// Q-learning
				double R = 1;
				if (judgeLose()) { R = -100; init(); }
				else if (pillarsX <= -pillarsWidth) { R = 15; initPillars(); }
				q.feedback(state, action, R);
			}flag = !flag;
			// draw
			if (time % 100000000 == 0 && q.greedy > 0)q.greedy -= 0.1;
			if (time++ > 1e9) {
				q.greedy = 0;
				drawPillars(pillarsX, windowHeight - pillarsY, &g);
				drawBird(birdX, windowHeight - birdY, &g);
				Delay(40);
				g.PicWrite("D:/LIGU.ppm");
				g.clear(0);
			}
		}
	}
	bool judgeLose() {
		if (birdY <= 0 || birdY >= windowHeight)return true;
		if (pillarsX < birdX && pillarsX + pillarsWidth > birdX && (birdY < pillarsY || birdY > pillarsY + holeHeight))
			return true;
		return false;
	}
	void drawBird(int x, int y, Graphics* g) {
		g->PaintSize = 10; g->PaintColor = 0xFFFF00;
		g->drawPoint(x, y);
		g->PaintSize = 0; g->PaintColor = 0xFFFFFF;
	}
	void drawPillars(int x, int y, Graphics* g) {
		g->PaintColor = 0x00FF00;
		g->drawRectangle(x, 0, x + pillarsWidth, y - holeHeight);
		g->drawRectangle(x, y, x + pillarsWidth, windowHeight);
		g->PaintColor = 0xFFFFFF;
	}
};

int main() {
	FlappyBird fb;
	fb.play();
}
