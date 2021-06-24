#ifndef _MATRIX_H
#define _MATRIX_H

#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<sstream>
#include<algorithm>
#include<iomanip>

using namespace std;

//自定义取绝对值的函数
double myabs(double x) { return x >= 0 ? x : -x; }

class Matrix
{
private:
	vector<vector<double>> mat;
	int rowSize;
	int colSize;
	

public:

	Matrix(int i, int j) :rowSize(i), colSize(j)
	{
		mat = vector<vector<double>>(i, vector<double>(j));
	}

	Matrix(vector<vector<double>>& matrix){
		mat = matrix;
		rowSize = matrix.size();
		colSize = matrix[0].size();
	}
	//构造函数，读入矩阵，并保存行数和列数
	Matrix(istream &in);

	//求转置矩阵
	vector<vector<double>> ReverseMat();

	//打印矩阵
	void matPrint();

	vector<double> fitParam();

	//获取行数和列数
	int getRowNum() const { return rowSize; }
	int getColNum() const { return colSize; }


	//行列交换
	void changeRow(int ri, int rj);
	void changeCol(int ci, int cj);

	//高斯全主元消元法
	void gaussEliminate();

	//LU 分解
	void LU(vector<vector<double>> &l, vector<vector<double>> &u);

	//求行列式的值
	double gaussDet();

	void inv(vector<vector<double>> &INV);


	//获取矩阵
	vector<vector<double>>& getMat()
	{
		return mat;
	}


	double getValue(int i, int j)
	{
		return mat[i - 1][j - 1];
	}


	//获取行元素
	vector<double> getRow(int i);

	//获取列元素
	vector<double> getCol(int j);


	//求两向量内积
	double dot(const vector<double> &v1, const vector<double> &v2)
	{
		double ret = 0;
		for (int i = 0; i < v1.size(); i++)
			ret += v1[i] * v2[i];

		return ret;
	}


	//计算向量长度

	double VecLen(const vector<double> &vec)
	{
		double ret=0;
		for (auto v : vec)
		{
			ret += v * v;
		}

		return sqrt(ret);
	}

	void Schmidt(vector<vector<double>> &q);

	void QR_Schmidt(vector<vector<double>> &Q, vector<vector<double>> &R);


	//矩阵相乘
	Matrix matMulti(const Matrix &mat2);

};

//高斯列主元求逆
void Matrix::inv(vector<vector<double>> &INV)
{
	int index = 0;
	for (auto &r : INV)
		r[index++] = 1;

	int i, j, k, rs, cs;
	double tmp, d;

	for (k = 0; k <= rowSize - 1; ++k)
	{
		d = 0.0;
		//选择该列最大元素
		for (i = k; i <= rowSize - 1; ++i)
		{
			if (fabs(mat[i][k]) > d)
			{
				d = mat[i][k];
				rs = i;
			}
		}

		//交换单位矩阵I对应的行和，并按照主元素归一化,注意下标+1对应

		if (k != rs)
		{
			changeRow(k, rs);
			std::swap(INV[k], INV[rs]);
		}

		tmp = mat[k][k];
		//消元
		for (auto &v : mat[k])
			v = v / tmp;
		for (auto &v : INV[k])
			v = v / tmp;

		for (i = 0; i <= rowSize - 1; i++)
		{
			if (i != k)
			{
				tmp = mat[i][k] / mat[k][k];
				for (j = 0; j <= colSize - 1; ++j)
				{
					mat[i][j] = mat[i][j] - tmp * mat[k][j];
					if (fabs(mat[i][j]) <= 1e-10)
						mat[i][j] = 0;

					INV[i][j] = INV[i][j] - tmp * INV[k][j];
				}
			}
		}
	}
	/*
	   for(i=0;i<=rowSize-2;i++)
	   {
		   tmp=mat[i][colSize-1]/mat[rowSize-1][colSize-1];
		   for(j=0;j<=colSize-1;j++)
		   {
			   INV[i][j]=INV[i][j]-tmp*INV[rowSize-1][j];
		   }
	   }
	*/
}

//高斯全主元化为上三角后求行列式的值
double Matrix::gaussDet()
{
	gaussEliminate();
	int ret = 1;
	int index = 0;
	for (auto v : mat)
	{
		ret *= v[index++];
	}
	return ret;
}

Matrix::Matrix(istream &in)
{
	string line;
	double word;
	vector<double> vec;
	while (getline(in, line))
	{
		istringstream record(line);
		while (record >> word)
			vec.push_back(word);

		mat.push_back(vec);
		colSize = vec.size();
		vec.clear();
	}
	rowSize = mat.size();
}

void Matrix::changeRow(int ri, int rj)
{
	std::swap(mat[ri], mat[rj]);
}

void Matrix::changeCol(int ci, int cj)
{
	for (auto &v : mat)
		std::swap(v[ci], v[cj]);
}

//高斯全主元消元法求上三角
void Matrix::gaussEliminate()
{
	int i, j, k, rs, cs;
	double tmp, d;

	for (k = 0; k <= rowSize - 2; ++k)
	{
		d = 0.0;
		for (i = k; i <= rowSize - 1; ++i)
		{
			for (j = k; j <= colSize - 1; ++j)
			{
				if (fabs(mat[i][j]) > d)
				{
					d = mat[i][j];
					rs = i;
					cs = j;
				}
			}
		}

		if (k != rs)
			changeRow(k, rs);
		if (k != cs)
			changeCol(k, cs);

		d = mat[k][k];

		for (i = k + 1; i <= rowSize - 1; ++i)
		{
			double tmp = mat[i][k] / mat[k][k];
			for (j = k; j <= colSize - 1; ++j)
			{
				mat[i][j] = mat[i][j] - tmp * mat[k][j];
				if (fabs(mat[i][j]) <= 1e-10)
					mat[i][j] = 0;
			}
		}
	}
}

/*
 * LU decomposition
 * LU(vector<vector<double>> &l,vector<vector<double>> &u)
 */

void Matrix::LU(vector<vector<double>> &l, vector<vector<double>> &u)
{
	int index = 0;
	for (auto &r : l)
	{
		r[index++] = 1;
	}

	for (auto v : mat)
		u.push_back(v);

	double d = 0.0;
	for (int k = 0; k < rowSize - 1; ++k)
	{
		for (int i = k + 1; i <= rowSize - 1; ++i)
		{
			l[i][k] = u[i][k] / u[k][k];
			for (int j = k; j <= colSize - 1; ++j)
			{
				u[i][j] = u[i][j] - l[i][k] * u[k][j];
				if (fabs(u[i][j]) < 1e-10)
					u[i][j] = 0;

			}

		}
	}
}


void Matrix::Schmidt(vector<vector<double>> &q)
{

	for (int j = 2; j <= colSize; j++)
	{
		double ret;
		double retk;
		for (int k = 1; k < j; k++)
		{
			ret = dot(getCol(j), getCol(k));
			retk = dot(getCol(k), getCol(k));

			for (auto &r : mat)
			{
				r[j - 1] = r[j - 1] - ret / retk * r[k - 1];
			}
		}
	}

	for (int j = 1; j <= colSize; j++)
	{
		double b = VecLen(getCol(j));
		for (auto &r : mat)
		{
			r[j - 1] = r[j - 1] / b;
		}
	}
	q = mat;
}


void Matrix::QR_Schmidt(vector<vector<double>> &Q, vector<vector<double>> &R)
{
	int index = 0;
	for (auto &vec : R)
	{
		vec[index++] = 1;
	}

	for (int j = 2; j <= colSize; j++)
	{
		double ret;
		double retk;
		for (int k = 1; k < j; k++)
		{
			ret = dot(getCol(j), getCol(k));
			retk = dot(getCol(k), getCol(k));
			R[k - 1][j - 1] = ret / retk;
			for (auto &r : mat)
			{
				r[j - 1] = r[j - 1] - ret / retk * r[k - 1];
			}
		}
	}

	for (int j = 1; j <= colSize; j++)
	{
		double b = VecLen(getCol(j));

		for (auto &v : R[j - 1])
		{
			v = v * b;
		}
		for (auto &r : mat)
		{
			r[j - 1] = r[j - 1] / b;
		}
	}
	Q = mat;

}

/**
 * vector<double> getRow(int i)
 * return the i row
 */

vector<double> Matrix::getRow(int i)
{
	return mat[i - 1];
}

/**
 * vector<double> getCol(int j)
 */
vector<double> Matrix::getCol(int j)
{
	vector<double> col;
	for (auto c : mat)
		col.push_back(c[j - 1]);
	return col;
}

/**
 * Matrix multiplication
 * matrix1=m*n   matrix2=n*k   matrix3=m*k
 */

Matrix Matrix::matMulti(const Matrix &matrix2)
{
	if (getColNum() != matrix2.getRowNum())
	{
		cout << "The col of mat1 is not equal the row of the mat2" << endl;
	}
	int m = getRowNum();
	int n = getColNum();
	int k = matrix2.getColNum();

	Matrix matrix3(m, k);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < k; ++j)
			for (int l = 0; l < n; ++l)
				matrix3.mat[i][j] += mat[i][l] * matrix2.mat[l][j];
	return matrix3;
}

/**
 * print the mat;
 */
void Matrix::matPrint()
{
	for (auto c : mat)
	{
		for (auto w : c)
		{
			cout << setiosflags(ios::fixed);
			cout << setw(10) << setprecision(7) << w;
		}
		cout << endl;
	}
}


//求转置矩阵
vector<vector<double>> Matrix::ReverseMat() {
	vector<vector<double>> Remat(mat[0].size(), vector<double>());
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			Remat[j].push_back(mat[i][j]);
		}
	}
	return Remat;
}


vector<double> Matrix::fitParam() {
	vector<double> C(6, 0);
	vector<double> res(6, 0);
	for (int i = 0; i < mat.size(); i++)
		C[i] = mat[i][0];
	double k1= (C[3]*C[3] - 4 * C[4]*C[5])*(C[3]*C[3]+C[4]*C[4]-2*C[4]*C[5]+C[5]*C[5]);
	double k2 = sqrt(1 - (C[3]*C[3] / (C[3]*C[3] + C[4]*C[4] - 2 * C[4]*C[5] + C[5]*C[5])));
	double k3= C[3]*C[3] - 4 * C[4]*C[5];

	int absC4 = (myabs(C[3]) / C[3]);
	int absC5_C6 = (myabs(C[4]-C[5]) / (C[4]-C[5]));

	double A0 = exp(0.5*(2*C[0] + 2 * (C[2]*C[2] * C[4] + C[1]*C[1]* C[5]) / k3 + (-2 * C[1]*C[2]*C[3]*C[3]*C[3] + C[1]*C[1]* C[3]*C[3] * C[4] - C[2]*C[2] * C[3]*C[3] * C[4] - C[1]*C[1] * C[3]*C[3] * C[5] + C[2]*C[2] * C[3]*C[3] * C[5]) / k1 + ((C[1]*C[1] - C[2]*C[2])*C[3] + 2 * C[1]*C[2]*(C[4] - C[5]))*k2*sqrt(1- k2*k2) / k3));
	
	double X0= (C[1]*C[3]*C[3] * (C[4] - C[5]) - C[2]*C[3]*C[3]*C[3]) / k1 + 2 * C[1]*C[5] / k3 + (C[2]*C[4]*k2 - C[1]*C[3]*k2 - C[2]*C[4]*k2)*sqrt(1 - k2*k2) / k3 * absC5_C6 * absC4;
	
	double Y0 = (C[2] * C[3] * C[3] * (C[5] - C[4]) - C[2] * C[3] * C[3] * C[3]) / k1 + 2 * C[2] * C[4] / k3 + (C[1] * C[5] * k2 + C[2] * C[3] * k2 - C[1] * C[4] * k2)*sqrt(1 - k2 * k2) / k3* absC5_C6 * absC4;
	
	double Delta_x= sqrt((C[4] + C[5] - myabs(C[4] - C[5])*k2 - myabs(C[3]) * sqrt(1 - k2 * k2)) / k3);

	double Delta_y= sqrt((C[4] + C[5] + myabs(C[4] - C[5])*k2 + myabs(C[3]) * sqrt(1 - k2 * k2)) / k3);

	double cos = absC4*(0.5*sqrt(2+2 * absC5_C6 * k2));
	res[0] = A0;
	res[1] = X0;
	res[2] = Y0;
	res[3] = Delta_x;
	res[4] = Delta_y;
	res[5] = cos;
	return res;
}

#endif


