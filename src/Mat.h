#ifndef _MAT_H
#define _MAT_H
template<class T>
class Mat
{
public:
	/* ---------------- 核心数据 ---------------- */
	T* data = NULL;
	int rows = 0, cols = 0;
	/* ---------------- 析构函数 ---------------- */
	~Mat(){ free(data); }
	/* ---------------- 零元 ---------------- */
	void zero(const int _rows, const int _cols) {
		free(data);
		data = (T*)malloc(sizeof(T) * _rows * _cols);
		memset(data, 0, sizeof(T) * _rows * _cols);
		rows = _rows;	cols = _cols;
	}
	/* ---------------- 单位元 ---------------- */
	void E(const int _rows) {
		zero(_rows, _rows);
		for (int i = 0; i < rows; i++) {
			*(data + i * cols + i) = 1;
		}
	}
	/* ---------------- 取值 赋值 ---------------- */
	void setValue(int i, int j, T value) {
		*(data + i * cols + j) = value;
	}
	T getValue(int i, int j) const {
		return *(data + i * cols + j);
	}
	/* ---------------- max min ---------------- */
	T max() const{
		T max = *data ;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				max = max >= *(data + i * cols + j) ? max : *(data + i * cols + j);
			}
		}
		return max;
	}
	T min() const{
		T min = getValue(0, 0);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				min = min <=  *(data + i * cols + j) ? min :  *(data + i * cols + j);
			}
		}
		return min;
	}
	/*----------------赋值 [ = ]----------------*/
	/*----------------加法 [ + ]----------------*/
	/*----------------乘法 [ * ]----------------*/
	void mult(const Mat& a, const Mat& b, Mat& ans) const {
		if (a.cols != b.rows) {
			return;
		}
		Mat ansTemp;
		ansTemp.zero(a.rows, b.cols);
		for (int i = 0; i < a.rows; i++) {
			for (int j = 0; j < b.cols; j++) {
				T sum;
				memset(&sum, 0, sizeof(sum));
				for (int k = 0; k < a.cols; k++) {
					T aV = *(a.data + i * a.cols + k);
					T bV = *(b.data + k * b.cols + j);
					sum += aV * bV;
				}
				ansTemp.setValue(i, j, sum);
			}
		}
		free(ans.data);
		ans.data = ansTemp.data; ansTemp.data = NULL;
	}
	/*----------------数乘 [ * ]----------------*/
	/*----------------转置 [ trans() ]----------------*/
	/*----------------余子式 [ comi ]----------------*/
	/*----------------取逆 [ inv ]----------------*/
	/*----------------行列式 [ abs() ]----------------*/
	/*--------------伴随矩阵 [ adj() ]----------------*/
	/*----------------特征值特征向量 [ eig() ]----------------*/
};
#endif