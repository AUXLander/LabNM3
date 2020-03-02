#include <iostream>
#include <algorithm>
#include <locale>
#include "TStructs.h"

using namespace std;

typedef long double (*Function)(long double x, long double y);

long double **Data;
Function ruleFunc[6];

Function funcDiff = [](long double x, long double y) {
	return 4.L;// -expl(-x * y * y);
};


//размерность сетки
unsigned int n;// = 4;
unsigned int m;// = 4;

unsigned int _n;// = n + 1;
unsigned int _m;// = m + 1;

unsigned int _nH;// = _n / 2;
unsigned int _mH;// = _m / 2;

long double X0 = 0.L;// 1.L;
long double Y0 = 0.L;// 2.L;
long double X1 = 2.L;
long double Y1 = 1.L;// 3.L;

long double  h;// = (X1 - X0) / n;
long double  k;// = (Y1 - Y0) / m;

unsigned int Smax = 10000;
long double accretion_min = 1e-008L;


long double Xi(unsigned int i) {
	return X0 + i * h;
}
long double Yj(unsigned int j) {
	return Y0 + j * k;
}

void InitRule(Function *functionArray) {
	//y change
	functionArray[0] = [](long double x, long double y) {
		return -powl(y - .5L, 2.L);// (y - 2)* (y - 3);
	};

	//x change
	functionArray[1] = [](long double x, long double y) {
		return .75L - powl(x - 1.L, 2.L);// x* (x - 1)* (x - 2);
	};

	//y change
	functionArray[2] = [](long double x, long double y) {
		return 1.L - powl(x - 1.L, 2.L) - powl(y - .5L, 2.L);
	};

	//x change
	functionArray[3] = [](long double x, long double y) {
		return 1.L - powl(x - 1.L, 2.L) - powl(y - .5L, 2.L);
	};

	//y change
	functionArray[4] = [](long double x, long double y) {
		return -powl(y - .5L, 2.L);// y* (y - 2)* (y - 3);
	};

	//x change
	functionArray[5] = [](long double x, long double y) {
		return .75L - powl(x - 1.L, 2.L);// (x - 1)* (x - 2);
	};
}

void InitGrid() {

	Data = new long double*[_n];
	for (int i = 0; i < _n; i++) {
		Data[i] = new long double[_m];

		for (int j = 0; j < _m; j++) {
			Data[i][j] = 0;
		}
	}

	for (int j = 0; j < _m; j++) {
		Data[0][j] = ruleFunc[0](Xi(0), Yj(j));
	}
	for (int i = 0; i < _nH; i++) {
		Data[i][m] = ruleFunc[1](Xi(i), Yj(m));
	}
	for (int j = _mH; j < _m ; j++) {
		Data[_nH][j] = ruleFunc[2](Xi(_nH), Yj(j));
	}
	for (int i = _nH; i < _n; i++) {
		Data[i][_mH] = ruleFunc[3](Xi(i), Yj(_mH));
	}
	for (int j = 0; j < _mH; j++) {
		Data[n][j] = ruleFunc[4](Xi(n), Yj(j));
	}
	for (int i = 0; i < _n; i++) {
		Data[i][0] = ruleFunc[5](Xi(i), Yj(0));
	}
}

void PrintGrid() {
	for (int j = m; j >= 0; j--) {
		printf("%d: ", j);

		for (int i = 0; i <= n; i++) {
			if (Data[i][j] >= 0) {
				if (Data[i][j] == 0) {
					if ((j <= m && j > _mH) && (i <= n && i >= _nH)) {
						printf("        \t");
					}
					else {
						printf(" %.6f\t", 0);
					}
				}
				else {
					printf(" %.6f\t", Data[i][j]);
				}
			}
			else {
				printf("%.6f \t", Data[i][j]);
			}
		}
		cerr << endl;
	}
	cerr << endl;
}

long double w = 1; // параметр метода (в интервале от 0 до 2)
long double h2;
long double k2;
long double a2;


// Calculate v[i][j]
long double Vspp(unsigned int i, unsigned int j, long double* accretion = nullptr) {
	long double Vspp = 0;
	long double Vs = Data[i][j];

	Vspp -= h2 * (Data[i + 1][j] + Data[i - 1][j]);
	Vspp -= k2 * (Data[i][j + 1] + Data[i][j - 1]);
	Vspp += a2 * Data[i][j] * ((1.L / w) - 1);
	Vspp += funcDiff(Xi(i), Yj(j));
	Vspp /= a2;
	Vspp *= w;

	if (accretion != nullptr) {
		(*accretion) = abs(Vspp - Vs);
	}

	return Vspp;
}

void TopRelaxationMethod(bool output = false) {
	h2 = -((long double)n / (long double)(X1 - X0)) * ((long double)n / (long double)(X1 - X0));
	k2 = -((long double)m / (long double)(Y1 - Y0)) * ((long double)m / (long double)(Y1 - Y0));
	a2 = -2 * (h2 + k2);

	unsigned int Scur = 0;
	long double accretion = 0;
	long double accretion_max = 0;

	do {
		accretion_max = 0;
		
		for (int j = _m - 2; j >= _mH; j--) {
			for (int i = 1; i < _nH; i++) {
				Data[i][j] = Vspp(i, j, &accretion);
				accretion_max = max(accretion, accretion_max);
			}
		}

		for (int j = _mH - 1; j > 0; j--) {
			for (int i = 1; i < n; i++) {
				Data[i][j] = Vspp(i, j, &accretion);
				accretion_max = max(accretion, accretion_max);
			}
		}

		Scur++;
	} while ((accretion_max >= accretion_min) && (Scur < Smax));

	if (output) {
		if (accretion_max < accretion_min) {
			cerr << "Выход по точности" << endl;
			cerr << "Совершено итераций: " << Scur << endl;
			cerr << "Достигнутая точность: " << accretion_max << endl;
			cerr << "Установленная точность: " << accretion_min << endl;
		}
		else {
			cerr << "Выход по итерации S: " << Scur << " > Smax = " << Smax << endl;
			cerr << "Достигнутая точность: " << accretion_max << endl;
		}
		cerr << endl;
	}
}



int main() {
	setlocale(LC_ALL, "Russian");

	cerr << "u(x,y) = 1 - (x - 1)^2 - (y - 0.5)^2;" << endl;

	cerr << "Введите размерность сетки Ox: ";	cin >> n;
	cerr << "Введите размерность сетки Oy: ";	cin >> m;
	cerr << "Введите точность метода: ";		cin >> accretion_min;
	cerr << "Введите ограничение шагов: ";		cin >> Smax;
	cerr << "Введите параметр w метода ВР: ";	cin >> w;

	_n = n + 1;
	_m = m + 1;

	_nH = _n / 2;
	_mH = _m / 2;

	h = (X1 - X0) / n;
	k = (Y1 - Y0) / m;

	InitRule(ruleFunc);
	InitGrid();

	cerr << "Исходная сетка:" << endl << endl;
	PrintGrid();

	TopRelaxationMethod(true);
	PrintGrid();

	cout << endl;
	system("pause");
}