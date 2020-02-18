#include <iostream>
#include "TStructs.h"
#include "Region.h"

typedef long double (*Function)(long double x, long double y);

long double **Data;
Point2D<long double>* Region;
Function ruleFunc[6];
Function funcDiff = [](long double x, long double y) {
	return 8.L;
};


//размерность сетки
unsigned int n = 4;
unsigned int m = 4;
long double  h = 1.L / n;
long double  k = 1.L / m;
long double X0 = 0;
long double Y0 = 0;


using namespace std;

long double Xi(unsigned int i) {
	return X0 + i * h;
}
long double Yj(unsigned int j) {
	return Y0 + j * h;
}

void InitRule(Point2D<long double>* region, Function *functionArray) {
	region = new Point2D<long double>[6];
	region[0] = { 0,0 };
	region[1] = { 0,2 };
	region[2] = { 1,2 };
	region[3] = { 1,1 };
	region[4] = { 2,1 };
	region[5] = { 2,0 };

	//y change
	functionArray[0] = [](long double x, long double y) {
		return 1.L;
	};
	//x change
	functionArray[1] = [](long double x, long double y) {
		return 2.L;
	};
	//y change
	functionArray[2] = [](long double x, long double y) {
		return 3.L;
	};
	//x change
	functionArray[3] = [](long double x, long double y) {
		return 4.L;
	};
	//y change
	functionArray[4] = [](long double x, long double y) {
		return 5.L;
	};
	//x change
	functionArray[5] = [](long double x, long double y) {
		return 6.L;
	};
}

void InitGrid(long double*** data, unsigned int _n, unsigned int _m) {

	_n = _n + 1;
	_m = _m + 1;
	
	(*data) = new long double*[_n];
	for (int i = 0; i < _n; i++) {
		(*data)[i] = new long double[_m];

		for (int j = 0; j < _m; j++) {
			(*data)[i][j] = 0;
		}
	}
	for (int j = 0; j < _m; j++) {
		(*data)[0][j] = ruleFunc[0](Xi(0), Yj(j));
	}
	for (int i = 0; i < n / 2; i++) {
		(*data)[i][m] = ruleFunc[1](Xi(i), Yj(m));
	}
	for (int j = m / 2; j < _m ; j++) {
		(*data)[n / 2][j] = ruleFunc[2](Xi(n / 2), Yj(j));
	}
	for (int i = n / 2; i < _n; i++) {
		(*data)[i][m / 2] = ruleFunc[3](Xi(i), Yj(m / 2));
	}
	for (int j = 0; j < m / 2; j++) {
		(*data)[n][j] = ruleFunc[4](Xi(n), Yj(j));
	}
	for (int i = 0; i < _n; i++) {
		(*data)[i][0] = ruleFunc[5](Xi(i), Yj(0));
	}
}


void TopRelaxationMethod(long double ***data) {
	
	long double w = 1;			// параметр метода (в интервале от 0 до 2)
	int Nmax = 10000;		// максимальное число итераций (не менее 1)
	int S = 0;				// счетчик итераций
	long double eps = 0.001; // минимально допустимый прирост
	long double eps_max = 0;		// текущее значение прироста
	long double eps_cur = 0;		// для подсчета текущего значения прироста
	long double a2, k2, h2;		// ненулевые элементы матрицы (–A)

	//double f[n + 1][m + 1]; // f( x,y) из дифф. уравнения в узлах сетки
	long double a = 0, b = 2, c = 0, d = 2;		// границы области


	long double v_old;		// старое значение преобразуемой компоненты вектора v
	long double v_new;		// новое значение преобразуемой компоненты вектора v
	bool f = false;			// условие остановки
	h2 = -(n / (b - a)) * (n / (b - a));
	k2 = -(m / (d - c)) * (m / (d - c));
	a2 = -2 * (1.L /(h*h)  + 1.L/(k*k));
	while (!f) {
		eps_max = 0;
		for (int j = 1; j < m; j++) {
			for (int i = 1; i < n; i++) {
				v_old = (*data)[i][j];
				v_new = -w * (h2 * ((*data)[i + 1][j] + (*data)[i - 1][j]) + k2 * ((*data)[i][j + 1] + (*data)[i][j - 1]));
				v_new = v_new + (1 - w) * a2 * (*data)[i][j] + w;// *f[i][j];
				v_new = v_new / a2;
				
				eps_cur = abs(v_old - v_new);

				if (eps_cur > eps_max) {
					eps_max = eps_cur;
				}

				(*data)[i][j] = v_new;
			}
		}
		S = S + 1;
		if ((eps_max < eps) || (S >= Nmax)) {
			f = true;
		}
	}
}



int main() {
	InitRule(Region, ruleFunc);
	InitGrid(&Data, n, m);
	TopRelaxationMethod(&Data);
	for (int j = 0; j <= m; j++) {
		for (int i = 0; i <= n; i++) {
			cout << Data[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
	cout << "Hello" << ruleFunc[0](0,0) << Data[0][0];
}