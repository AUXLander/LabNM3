#pragma once
#include <math.h>

template<class T1>
struct Point2D {
	T1 min;
	T1 max;
};

template<class T1>
struct Square2D {
	Point2D<T1> X;
	Point2D<T1> Y;
};