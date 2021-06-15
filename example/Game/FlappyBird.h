#ifndef  FLAPPY_BIRD_H
#define  FLAPPY_BIRD_H
/******************************************************************************
*                    FlappyBird
******************************************************************************/
class FlappyBird {
public:
	int W, H;
	FlappyBird(int _W = 400, int _H = 600) { W = _W, H = _H; init();}
	int birdY = H / 2, birdX = 100, pillarsWidth = 80, holeHeight = 100, pillarsX = W, pillarsY = 0;
	int pillarsSpeed = 10, birdSpeed = 20;
	void init() {
		birdY = H / 2; initPillars();
		pillarsX = W;
		pillarsY = rand() / double(RAND_MAX) * (H - 100 - 100) + 100;	//[st,ed)
	}
	void play(bool action) {
		pillarsX -= pillarsSpeed;
		birdY    -= birdSpeed;
		if (action) birdY += 3 * birdSpeed;
		//judgeLose
		if (birdY <= 0 || birdY >= H
		|| (pillarsX < birdX && pillarsX + pillarsWidth > birdX && (birdY < pillarsY || birdY > pillarsY + holeHeight))
		) return true;
		return false;
	}
};
#endif