#ifndef _MAT_H
#define _MAT_H
template<class T>
class Mat
{
public:
	T* data = NULL;
	int rows = 0, cols = 0;
	Mat() { ; }
	Mat(const int _rows) {
		E(_rows);
	}
	Mat(const int _rows, const int _cols) {
		zero(_rows, _cols);
	}
	~Mat(){ free(data); }
	/* ---------------- 零元 ---------------- */
	void zero(const int _rows, const int _cols) {
		free(data);
		data = (T*)malloc(sizeof(T) * _rows * _cols);
		memset(data, 0, sizeof(T) * _rows * _cols);
		rows = _rows;
		cols = _cols;
	}
	/* ---------------- 单位元 ---------------- */
	void E(const int _rows) {
		zero(_rows, _rows);
		for (int i = 0; i < rows; i++) {
			setValue(i, i, 1);
		}
	}
	/* ---------------- 取值 赋值 ---------------- */
	void setValue(int i, int j, T value) {
		*(data + i * cols + j) = value;
	}

	T getValue(int i, int j) const {
		return *(data + i * cols + j);
	}
	/*----------------赋值 [ = ]----------------*/
	void operator= (const Mat& b) {
		zero(b.rows, b.cols);
		memcpy(data, b.data, sizeof(T) * b.rows * b.cols);
		//return *this会出现乱码???
	}
	/*----------------加法 [ + ]----------------*/
	Mat operator+ (const Mat& b) const {
		if (rows != b.rows || cols != b.cols) {
			exit(-1);
		}
		Mat Temp(rows, b.cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				T aV = getValue(i, j);
				T bV = b.getValue(i, j);
				Temp.setValue(i, j, aV + bV);
			}
		}
		return Temp;
	}
	Mat operator+ (const T b) const {
		Mat Temp(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				T aV = getValue(i, j);
				Temp.setValue(i, j, aV + b);
			}
		}
		return Temp;
	}
	friend Mat operator+ (const T b, const Mat a) {
		return a + b;
	}
	/*----------------乘法 [ * ]----------------*/
	void mult(const Mat& a, const Mat& b, Mat& ans) const {
		if (a.cols != b.rows) {
			exit(-1);
		}
		Mat ansTemp(a.rows, b.cols);
		for (int i = 0; i < a.rows; i++) {
			for (int j = 0; j < b.cols; j++) {
				T sum;
				memset(&sum, 0, sizeof(sum));
				for (int k = 0; k < a.cols; k++) {
					T aV = a.getValue(i, k);
					T bV = b.getValue(k, j);
					sum += aV * bV;
				}
				ansTemp.setValue(i, j, sum);
			}
		}
		free(ans.data);
		ans.data = ansTemp.data; ansTemp.data = NULL;
	}
	Mat operator* (const Mat& b) const {
		if (cols != b.rows) {
			exit(-1);
		}
		Mat Temp(rows, b.cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < b.cols; j++) {
				T sum;
				memset(&sum, 0, sizeof(sum));
				for (int k = 0; k < cols; k++) {
					T aV = getValue(i, k);
					T bV = b.getValue(k, j);
					sum += aV * bV;
				}
				Temp.setValue(i, j, sum);
			}
		}
		return Temp;
	}
	/*----------------数乘 [ * ]----------------*/
	Mat operator* (const T b) const {
		Mat Temp(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				T aV = getValue(i, j);
				Temp.setValue(i, j, aV * b);
			}
		}
		return Temp;
	}
	friend Mat operator* (const T b, const Mat a) {
		return a * b;
	}
	/*----------------转置 [ trans() ]----------------*/
	Mat trans() const {
		Mat Temp(cols, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				T t = getValue(i, j);
				Temp.setValue(j, i, t);
			}
		}
		return Temp;
	}
	/*----------------余子式 [ comi ]----------------*/
	Mat comi(const int x,const int y) const {
		Mat Temp(rows - 1, cols - 1);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (i == x || j == y)continue;
				T t = getValue(i, j);
				Temp.setValue(i, j, t);
			}
		}
		if ((x + y) % 2 == 1)
			Temp = -1 * Temp;
		return Temp;
	}
	/*----------------取逆 [ inv ]----------------*/
	Mat inv() const {
		Mat Temp;
		Temp = this->trans();
		return Temp;
	}
	/*----------------行列式 [ abs() ]----------------*/
	friend T abs(Mat a) {
		if (a.cols != a.rows) {
			exit(-1);
		}
		if (a.rows == 1) {
			return a.getValue(0, 0);
		}
		T sum;
		memset(&sum, 0, sizeof(sum));
		for (int i = 0; i < a.rows; i++) {
			printf("%d %d\n",i, a.rows);
			Mat t = a.comi(i, 0);
			sum += a.getValue(0, i) * abs(t);
		}
		return sum;
	}
	/*--------------伴随矩阵 [ adj() ]----------------*/
	Mat adj(Mat a) {

	}
	/*----------------特征值特征向量 [ eig() ]----------------*/
};
#endif