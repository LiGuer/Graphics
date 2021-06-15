#ifndef  SNAKE_H
#define  SNAKE_H
#include <conio.h>
#include <stdlib.h>
#include <stdio.h>
/******************************************************************************
*                    Snake 贪吃蛇
*	[算法]: 环形链表
-------------------------------------------------------------------------------
*	[Example]:
	Snake S(50, 50);
	while (S.play()) { 
		Snake::SnakeCell* cur = S.snake;
		while (cur->next != S.snake) { drawRectangle(cur->x, cur->y); cur = cur->next; }
		drawRectangle(cur->x, cur->y);
		drawRectangle(S.apple.x, S.apple.y);
		Delay(80);
	}
******************************************************************************/
class Snake {
public:
	struct SnakeCell { 
		int x = 0, y = 0; SnakeCell* next = NULL, * prev = NULL;
		void init(int _x, int _y, SnakeCell* _next, SnakeCell* _prev) { x = _x; y = _y; next = _next; prev = _prev;}
		void init(int _x, int _y) { x = _x; y = _y; }
	};
	int H, W, dH = 0, dW = 1;
	SnakeCell* snake, apple;
	Snake(int _H,int _W) {
		H = _H;
		W = _W;
		snake = new SnakeCell; snake->init(H / 2, W / 2, snake, snake);
		apple.init(rand() % H, rand() % W);
	}
	bool play() {
		interactive();
		if (snake->x < 0 || snake->x >= H 
		||  snake->y < 0 || snake->y >= W
		||  eatSelf()) return false;
		else if (snake->x == apple.x && snake->y == apple.y) {
			SnakeCell* head = new SnakeCell;
			head->init(
				snake->x + dH,
				snake->y + dW,
				snake,
				snake->prev
			);
			snake = snake->prev = snake->prev->next = head;
			apple.init(rand() % H, rand() % W);
		}
		else {
			snake->prev->init(
				snake->x + dH,
				snake->y + dW
			); snake = snake->prev;
		}
		return true;
	}
	bool eatSelf() {
		Snake::SnakeCell* cur = snake->next;
		while (cur != snake)
			if (cur->x == snake->x && cur->y == snake->y) return true;
			else cur = cur->next;
		return false;
	}
	static char Snake_AI(SnakeCell* snake, SnakeCell& apple, int H, int W) {
		// 10x10 可行
		static bool flag = 0;
		if (snake->x == H - 1 && snake->y == W - 2) return ' ';
		if (snake->x == H - 1 && snake->y == W - 1) return 's';
		if (snake->x == 0     && snake->y == W - 1) return 'a';
		if (snake->x == 0     && snake->y == W - 2) return ' ';
		if (flag == 0		  && snake->y == W - 2) { flag = 1; return 'w'; }
		if (flag == 1		  && snake->y == W - 2) { flag = 0; return 'a'; }
		if (flag == 0		  && snake->y == 0)     { flag = 1; return 'w'; }
		if (flag == 1		  && snake->y == 0)     { flag = 0; return 'd'; }
		return ' ';
	}
	void interactive() {
		int ch;
		if (_kbhit()) {
			ch = _getch();
			if (ch == 'd') dH = 0, dW = 1;
			if (ch == 'a') dH = 0, dW =-1;
			if (ch == 'w') dH = 1, dW = 0;
			if (ch == 's') dH =-1, dW = 0;
		}
		/*{
			ch = Snake_AI(snake, apple, H, W); printf("%d %d %c\n", snake->x, snake->y, ch);
			if (ch == 'd') dH = 0, dW = 1;
			if (ch == 'a') dH = 0, dW =-1;
			if (ch == 'w') dH = 1, dW = 0;
			if (ch == 's') dH =-1, dW = 0;
		}*/
		
	}
};
#endif // ! SNAKE_H
