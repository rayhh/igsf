#include<iostream>
#include<fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//mac下M_PI在math.h中有宏定义，所以这里我们选择行的宏定义
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIZE 1024
#define VALUE_MAX 1000

using namespace std;
////////////////////////////////////////////////////////////////////
//定义一个复数结构体
///////////////////////////////////////////////////////////////////
struct Complex_{
    double real;
    double imagin;
};
typedef struct Complex_ Complex_;

void Add_Complex(Complex_ * src1,Complex_ *src2,Complex_ *dst);
void Sub_Complex(Complex_ * src1,Complex_ *src2,Complex_ *dst);
void Multy_Complex(Complex_ * src1,Complex_ *src2,Complex_ *dst);
void getWN(double n,double size_n,Complex_ * dst);
void readSrcData(double(&arr)[12][12], double* xData, double* yData, double xShift, double yShift, int length, bool isBeam);
//void setInput(double* data, int size_n);
int FFT_remap(double * src,int size_n);
void FFT_func(double * src,Complex_ * dst,int size_n);
void fftAmp(Complex_ * dst, double* amp);
void gauss(double(&y)[SIZE], double sigma);
double findMaxFreq(double(&Freq)[SIZE]);
