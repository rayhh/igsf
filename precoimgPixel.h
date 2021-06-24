#pragma once
#ifndef _GAUSS_FIT_

#define _GAUSS_FIT_


#include<iostream>
#include<fstream>
#include<string>
#include<opencv2/highgui.hpp>
#include<opencv2/core.hpp>
#include<opencv2/imgcodecs.hpp>
#include<vector>
#include<cmath>
#include"cadef.h"

#define PIX_SIZE 3.8
#define PARTIClE_NUM 20000


using namespace cv;
using namespace std;





class Preco_imgPixel {
	private:
		Mat imagePixel;
		double weight;
	public:
		Preco_imgPixel(int row, int col, Mat& image,double w);
		~Preco_imgPixel() {}


		void showImage();
		vector<int> find_Max();
	    vector<vector<double>> CreateCoorFit(short max, bool isBeam);
		vector<vector<double>> CreateA(vector<double>& Fx);
	};


void coor_Convert(const char* path, double(&arr)[12][12],double theta);
void output_to_File(const char* outpath, double(&arr)[12][12]);
int channel_to_Pixel(chid pvID, double (&pChanel)[135], double(&arr)[12][12]);
vector<chid> epics_init();
double find_arr_max_double(double (&arr)[12][12]);
void Pixel_to_GrayScale(double maxChanel, double (&arr)[12][12]);
void output_fitvalue(const char* outpath, vector<double> &fitValue);

double AverageRandom(double min, double max);
double Normal2(double x, double y, double miu1, double miu2, Mat matrix, float detCov);
Point2i NormalRandom2(double miu1, double miu2, Mat matrix, float detCov, double min, double max);
void covMatrix(float a, float b, float c, float d, Mat& inverse_matrix, float& detCov);
void coor_Convert_Sim_Version(double(&arr)[12][12],double theta);
#endif // !_GAUSS_FIT_
