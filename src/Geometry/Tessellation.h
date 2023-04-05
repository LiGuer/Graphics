#include <vector>
#include "../../../Math/Matrix/Matrix.h"


// Regular hexagon Tessellation
void HexagonTessellation(int numX, int numY, double gridLen, std::vector<Mat<>>& points, std::vector<Mat<>>& edges) {
    points.clear();
    edges. clear();

    double
        x0 = 0,
        y0 = 0;
    Mat<> p(2), e(2, 2);

    for (int i = 0; i < numY; i++) {
        if (i != 0) {
            y0 -= 3.0 / 2 * gridLen; 

            if (i % 2 == 0)
                x0 += sqrt(3) / 2 * gridLen;
            else 
                x0 -= sqrt(3) / 2 * gridLen;
        }

        double x = x0;

        for (int j = 0; j < numX; j++) {
            if (j != 0)
                x += sqrt(3) * gridLen;

            points.push_back(p = { x, y0 });
            edges. push_back(e = { x - sqrt(3) / 2 * gridLen, x - sqrt(3) / 2 * gridLen, y0 - gridLen / 2, y0 + gridLen / 2 });
            edges. push_back(e = { x - sqrt(3) / 2 * gridLen, x, y0 + gridLen / 2, y0 + gridLen });
            edges. push_back(e = { x, x + sqrt(3) / 2 * gridLen, y0 + gridLen, y0 + gridLen / 2 });
            
            if (j == 0) {
                edges.push_back(e = { x, x - sqrt(3) / 2 * gridLen, y0 - gridLen, y0 - gridLen / 2 });
            }
            if (j == numX - 1) {
                if (i % 2 == 0)
                    edges.push_back(e = { x + sqrt(3) / 2 * gridLen, x, y0 - gridLen / 2, y0 - gridLen });
                edges.push_back(e = { x + sqrt(3) / 2 * gridLen, x + sqrt(3) / 2 * gridLen, y0 + gridLen / 2, y0 - gridLen / 2 });
            }
            if (i == numY - 1) {
                edges.push_back(e = { x + sqrt(3) / 2 * gridLen, x, y0 - gridLen / 2, y0 - gridLen });
                edges.push_back(e = { x, x - sqrt(3) / 2 * gridLen, y0 - gridLen, y0 - gridLen / 2 });
            }
        }
    }
}

// Random Rectangle Tessellation
void RandomRectangleTessellation(
	Mat<>& area, std::vector<Mat<>>& rects, int N, int sizeThreshold = 0, int randThreshold = 0
) {
	Mat<> p(4), p1(4), p2(4);

	rects.clear();
	rects.push_back(area);

	for (int i = 0; i < N; i++) {
		int n = rects.size();

		for (int j = 0; j < n; j++) {
			p1 = p2 = rects[j];

			int ind = (i % 2 == 0) ? 0 : 1;

			if (rects[j][ind + 2] > sizeThreshold) {
				double ra = rand() / (double)RAND_MAX * (1 - 2 * randThreshold) + randThreshold;
				ra *= rects[j][ind + 2];

				p1[ind + 2] = ra;
				p2[ind + 2] = rects[j][ind + 2] - ra;
				p1[ind] = rects[j][ind] - rects[j][ind + 2] / 2 + ra / 2;
				p2[ind] = rects[j][ind] + ra / 2;

				rects[j] = p1;
				rects.push_back(p2);
			}
		}
	}
}

