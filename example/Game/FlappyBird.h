#ifndef  FLAPPY_BIRD_H
#define  FLAPPY_BIRD_H
/******************************************************************************
*                    FlappyBird
******************************************************************************/
class FlappyBird {
public:
	int W, H, birdX, birdY, pillarsX, pillarsY,
		pillarsW = 80, holeH = 100, pillarsSpeed = 10, birdSpeed = 20;
	FlappyBird(int _W = 400, int _H = 600) { W = _W, H = _H; init();}
	void init() {
		birdY = H / 2, birdX = 100;
		pillarsX = W;
		pillarsY = rand() / double(RAND_MAX) * (H - 100 - 100) + 100;
	}
	void play(bool action) {
		pillarsX -= pillarsSpeed;
		birdY    +=   -birdSpeed + (action ? 3 * birdSpeed : 0);
		if (pillarsX < 0) 
			pillarsX = W, 
			pillarsY = rand() / double(RAND_MAX) * (H - 100 - 100) + 100;
		//judgeLose
		if (birdY <= 0 || birdY >= H
		||((birdX > pillarsX && birdX < pillarsX + pillarsW) 
		&& (birdY < pillarsY || birdY > pillarsY + holeH))) return true;
		return false;
	}
};
#endif