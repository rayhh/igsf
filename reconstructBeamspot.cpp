#include<time.h>
#include<unistd.h>
#include<stdlib.h>
#include<vector>
#include<fstream>
#include<pthread.h>
#include<opencv2/imgproc.hpp>
#define CVUI_IMPLEMENTATION
#include "cvui.h"
#include "precoimgPixel.h"
#include "myMatrix.h"
#include "fft.h"


const double PI=3.1415926;           //定义PI 常量
bool myepicsStatus = false;             //定义与EPICS设备状态有关的变量
bool isBeam = false;                 //定义isBeam用于区别仿真调试模式和束流模式，默认为仿真调试模式
vector<chid> pvID;                   //定义过程变量通道ID
int imgshowPeriod = 102;             //显示更新周期
int igsfPeriod = 55*1000;                 //IGSF计算周期
int fftPeriod = 45*1000;                  //FFT计算周期
int epicsPeriod = 30*1000;                //EPICS通信周期



//******************************************加锁全局变量创建区域*******************************
//*******************************************************************************************
//*******************************************************************************************
//创建数组用于统计每个像素点粒子数目，共12*12=144个像素点
struct lockArr{
	pthread_mutex_t arr_lock;
	double arr[12][12];
};

//定义噪声门限值
struct lockNoiseThreshold{
	pthread_mutex_t noiseThreshold_lock;
	double nosieThreshold;
};

//定义互斥访问的EPICS相关变量，包括过程变量ID和过程
struct lockPV{
	pthread_mutex_t PV_lock;
	double pvChanel[135];
	bool px3Status;
	double hvSet;
	double hvRead;
	bool hvEnable;
	double integratedTimeSet;
	double integratedTimeRead;
	double px3Tempature;
	double px3Pressure;
	double px3Humidity;
};

//拟合结果存储区域
struct lockIGSFResult{
	pthread_mutex_t IGSFResult_lock;
	vector<double> fit_value;
};

//FFT计算结果存储区域
struct lockFFT{
	pthread_mutex_t FFT_lock;
	double fftOutputAmpX[SIZE];
	double fftOutputAmpY[SIZE];
	double fftOutputAmpMajor[SIZE];
	double fftOutputAmpMinor[SIZE];
	double fftMaxX;
	double fftMaxY;
	double fftMaxMajor;
	double fftMaxMinor; 
};

struct lockArr Arr = {PTHREAD_MUTEX_INITIALIZER, {0}};
struct lockNoiseThreshold NoiseThreshold = {PTHREAD_MUTEX_INITIALIZER, 0.05};
struct lockPV PV = {PTHREAD_MUTEX_INITIALIZER, {0}, false, 0, 0, false, 0, 0, 0, 0, 0};
struct lockIGSFResult IGSFResult = {PTHREAD_MUTEX_INITIALIZER,{0,0,0,0,0,0}};
struct lockFFT FFT = {PTHREAD_MUTEX_INITIALIZER, {0},{0},{0},{0},0,0,0,0};
//******************************************************************************************
//******************************************************************************************
//********************************互斥访问区域定义完毕****************************************



void* imshowWorker(void*) {
	//定义粒子源文件路径和像素文件存储的路径     
	const char* inpath = "../data/proton.txt";
	const char* outpath = "../data/pixelData.txt";
	//定义变量用于存储像素值中的最大值
	double maxChanel;
	//定义噪声门限值的上下限
	double noiseLow = 0.01;
	double noiseHigh = 0.2;
	//定义高压值的上下限
	double hvHigh = 1500.0;
	double hvLow = 0.0;
	//定义积分时间的上下限
	double integratedtimeLow = 1000.0;
	double integratedtimeHigh = 2000.0;
	//定义布尔变量maxConnDomain用于判决是否采用主连通域方法处理像素数据，默认为true
	bool maxConnDomain = true;
	//定义布尔变量用于判决是否勾勒拟合出的椭圆的轮廓，轮廓的大小默认为3*sigma，可在
	//(coefLow, coefHigh)范围内调节。
	bool isEllipse = false;
	double coefficient = 1.5;
	double coefLow = 1.0;
	double coefHigh = 3.0;
	//与显示有关的参数设置
	int pixelSize = 35;
	int beamspotStartX = 10;
	int beamspotStartY = 9;
	//定义Mat类画布矩阵并初始化尺寸
	Mat dst;
	dst.create(640,1010,CV_8UC3);
	dst = imread("/mnt/d/background.bmp");
	//定义线程内的局部变量imgArr  fitValue
	double imgArr[12][12] = {0};
	vector<double> fitValue(6,0);
	double noiseThreshold_temp = 0.05;
	//定义线程内关于EPICS的局部变量
	bool px3Status_temp;
	double hvRead_temp, integratedTimeRead_temp, px3Tempature_temp;
	double px3Pressure_temp, px3Humidity_temp;
	double integratedTimeSet_temp = 0;
	bool HVEnable_temp = 0;
	double HVSet_temp = 0;
	//定义与FFT有关的局部变量151
	double fftOutputAmpX_temp[SIZE], fftOutputAmpY_temp[SIZE], fftOutputAmpMajor_temp[SIZE];
	double fftOutputAmpMinor_temp[SIZE], fftMaxX_temp, fftMaxY_temp, fftMaxMajor_temp, fftMaxMinor_temp;
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////程序开始部分///////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//打开粒子源文件，并模拟像素电离室的粒子收集过程，构造出12*12的像素矩阵imgArr
	if(!isBeam){
		coor_Convert(inpath, imgArr, 0);
		//coor_Convert_Sim_Version(imgArr,0);
		maxChanel = find_arr_max_double(imgArr);
	}
	//画图的初始化以及界面初始化
	namedWindow("IGSF");
	cvui::init("IGSF");
	while(true){
		//将pv变量的通道数据映射成12*12的像素矩阵arr
		pthread_mutex_lock(&PV.PV_lock);
		px3Status_temp = PV.px3Status;
		px3Pressure_temp = PV.px3Pressure;
		px3Humidity_temp = PV.px3Humidity;
		px3Tempature_temp = PV.px3Tempature;
		hvRead_temp = PV.hvRead;
		integratedTimeRead_temp = PV.integratedTimeRead;
		PV.integratedTimeSet = integratedTimeSet_temp;
		PV.hvSet = HVSet_temp;
		PV.hvEnable = HVEnable_temp;
		pthread_mutex_unlock(&PV.PV_lock);
		if(isBeam){
			//myepicsStatus = channel_to_Pixel(pvID[0], pvChanel, arr);
			//output_to_File(outpath, arr);
			//像素矩阵arr的归一化处理以及负数噪声取反变成正值
			////////////test用途/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			ifstream in;
			in.open(outpath,ios::in);
			if(!in.is_open()){
				cout<<"open pixelData.txt failed"<<endl;
			}else{
				for(int i=0; i<12; i++){
					in>>imgArr[i][0]>>imgArr[i][1]>>imgArr[i][2]>>imgArr[i][3]>>imgArr[i][4]>>imgArr[i][5]>>imgArr[i][6]>>imgArr[i][7]>>imgArr[i][8]>>imgArr[i][9]>>imgArr[i][10]>>imgArr[i][11];
				}
				in.close();
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			maxChanel = find_arr_max_double(imgArr);
			Pixel_to_GrayScale(maxChanel, imgArr);
			maxChanel = find_arr_max_double(imgArr);
		}
		pthread_mutex_lock(&Arr.arr_lock);
		for(int i = 0; i<12; i++){
			for(int j=0; j<12; j++){
				Arr.arr[i][j] = imgArr[i][j];
			}
		}
		pthread_mutex_unlock(&Arr.arr_lock);
		pthread_mutex_lock(&IGSFResult.IGSFResult_lock);
		fitValue = IGSFResult.fit_value;
		pthread_mutex_unlock(&IGSFResult.IGSFResult_lock);
		pthread_mutex_lock(&NoiseThreshold.noiseThreshold_lock);
		NoiseThreshold.nosieThreshold = noiseThreshold_temp;
		pthread_mutex_unlock(&NoiseThreshold.noiseThreshold_lock);
		pthread_mutex_lock(&FFT.FFT_lock);
		for(int i = 0; i<SIZE; i++){
			fftOutputAmpMajor_temp[i] = FFT.fftOutputAmpMajor[i];
			fftOutputAmpMinor_temp[i] = FFT.fftOutputAmpMinor[i];
			fftOutputAmpX_temp[i] = FFT.fftOutputAmpX[i];
			fftOutputAmpY_temp[i] = FFT.fftOutputAmpY[i];
		}
		fftMaxX_temp = FFT.fftMaxX;
		fftMaxY_temp = FFT.fftMaxY;
		fftMaxMajor_temp = FFT.fftMaxMajor;
		fftMaxMinor_temp = FFT.fftMaxMinor;
		pthread_mutex_unlock(&FFT.FFT_lock);
		//////////////////////////////画图区域///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//填充dst的束斑部分
		dst = imread("/mnt/d/background.bmp");
		for (int i = 0; i < 12; i++){
			for (int j = 0; j < 12; j++){
				if(isBeam){
					cv::rectangle(dst,Point(beamspotStartX+j*pixelSize,beamspotStartY+i*pixelSize),Point(beamspotStartX+j*pixelSize+pixelSize,i*pixelSize+pixelSize+beamspotStartY),Scalar(imgArr[i][j]/maxChanel*255,imgArr[i][j]/maxChanel*255,imgArr[i][j]/maxChanel*255),-1);
				}
				else{
					cv::rectangle(dst,Point(beamspotStartX+i*pixelSize,beamspotStartY+pixelSize*12-j*pixelSize),Point(beamspotStartX+i*pixelSize+pixelSize,beamspotStartY-j*pixelSize-pixelSize),Scalar(imgArr[i][j]/maxChanel*255,imgArr[i][j]/maxChanel*255,imgArr[i][j]/maxChanel*255),-1);
				}
				
			}
		}
		cv::line(dst, Point(beamspotStartX,beamspotStartY+pixelSize*6), Point(beamspotStartX+pixelSize*12,beamspotStartY+pixelSize*6),Scalar(255,0,0));
		cv::line(dst, Point(beamspotStartX+pixelSize*6,beamspotStartY), Point(beamspotStartX+pixelSize*6,beamspotStartY+pixelSize*12),Scalar(255,0,0));
		cv::circle(dst, Point(fitValue[1]*pixelSize/3.8+pixelSize*6+beamspotStartX,-fitValue[2]*pixelSize/3.8+pixelSize*6+beamspotStartY),8,Scalar(255,0,0),-1);
		//拟合结果显示区
		cv::rectangle(dst, Point(110,442),Point(228,475),Scalar(255,255,255),-1);
		cv::rectangle(dst, Point(110,478),Point(228,509),Scalar(255,255,255),-1);
		cv::rectangle(dst, Point(110,512),Point(228,542),Scalar(255,255,255),-1);
		cv::rectangle(dst, Point(110,545),Point(228,577),Scalar(255,255,255),-1);
		cv::rectangle(dst, Point(110,580),Point(228,610),Scalar(255,255,255),-1);
		cvui::printf(dst,115,453,0.5,0xff0000,"%.3fmm",fitValue[3]);
		cvui::printf(dst,115,489,0.5,0xff0000,"%.3fmm",fitValue[4]);
		cvui::printf(dst,115,522,0.5,0xff0000,"%.3fmm",fitValue[1]);
		cvui::printf(dst,115,556,0.5,0xff0000,"%.3fmm",fitValue[2]);
		cvui::printf(dst,115,590,0.5,0xff0000,"%.2fDeg",acos(fitValue[5])*180/PI);
		//PX3设备参数设置区域
		if(myepicsStatus == true){
			cv::rectangle(dst, Point(868,361),Point(984,393),Scalar(0,255,0),-1);
		}
		else{
			cv::rectangle(dst, Point(868,361),Point(984,393),Scalar(0,0,255),-1);
		}
		if(px3Status_temp == true){
			cv::rectangle(dst, Point(868,396),Point(984,424),Scalar(0,255,0),-1);
		}
		else{
			cv::rectangle(dst, Point(868,396),Point(984,424),Scalar(0,0,255),-1);
		}
		cvui::printf(dst, 872, 468, 0.5, 0x000000, "%.2f",hvRead_temp);
		cvui::trackbar(dst, 445, 373, 290, &HVSet_temp, hvLow, hvHigh, 1, "%.1Lf", 0, 1.0);
		cvui::checkbox(dst, 459,428,"", &HVEnable_temp);		
		cvui::trackbar(dst, 445, 478, 290, &integratedTimeSet_temp, integratedtimeLow, integratedtimeHigh, 1, "%.1Lf", 0, 1.0);		
		cvui::printf(dst, 872, 500, 0.5, 0x000000, "%f", integratedTimeRead_temp);
		cvui::printf(dst, 872, 531, 0.5, 0x000000, "%f", px3Tempature_temp);
		cvui::printf(dst, 872, 563, 0.5, 0x000000, "%f", px3Pressure_temp);
		cvui::printf(dst, 872, 597, 0.5, 0x000000, "%f", px3Humidity_temp);
		//拟合参数设置区域
		cvui::checkbox(dst,243,440,"",&maxConnDomain);
		cvui::checkbox(dst,243,460,"",&isBeam);
		cvui::checkbox(dst,243,480,"",&isEllipse);
		cvui::printf(dst,243,505,0.4,0xffffff,"");
		cvui::trackbar(dst,243,520,180,&noiseThreshold_temp,noiseLow,noiseHigh,1, "%.2Lf",0, 0.01);
		cvui::printf(dst,243,575,0.4,0xffffff,"");
		cvui::trackbar(dst,243,590,180,&coefficient,coefLow,coefHigh,1, "%.2Lf",0, 0.1);
		
		if(isEllipse){
				ellipse(dst,Point(int(fitValue[1]*pixelSize/3.8+pixelSize*6+beamspotStartX),int(-fitValue[2]*pixelSize/3.8+pixelSize*6+beamspotStartY)), Size(int(fitValue[3]*coefficient*pixelSize/3.8), int(fitValue[4]*coefficient*pixelSize/3.8)),180-acos(fitValue[5])*180/PI,0,360,Scalar(0,0,255),1,8);
		}
		//傅里叶变换区域
		for(int i=0; i<9; i++){
			cvui::printf(dst,456,247-24-24*i,0.28,0x000000,"%.1f",(i+1)*0.1);
			cv::line(dst,Point(486,249-24-24*i),Point(490,249-24-24*i),Scalar(0,0,0));
		}
		for(int i=0; i<6; i++){
			cvui::printf(dst,486+69+69*i,251,0.28,0x000000,"%.2f",(i+1)*0.02);
			cv::line(dst,Point(486+78+78*i,249),Point(486+78+78*i,245),Scalar(0,0,0));
		}
		for(int i=0; i<512; i++){
			cv::circle(dst, Point(486+i,249-240*fftOutputAmpY_temp[i]),1,Scalar(255,0,0),-1);
			cv::circle(dst, Point(486+i,249-240*fftOutputAmpX_temp[i]),1,Scalar(0,0,255),-1);
			cv::circle(dst, Point(486+i,249-240*fftOutputAmpMajor_temp[i]),1,Scalar(0,0,0),-1);
			cv::circle(dst, Point(486+i,249-240*fftOutputAmpMinor_temp[i]),1,Scalar(128,255,128),-1);
		}
		cv::rectangle(dst, Point(443,307), Point(994,337), Scalar(0,0,0),-1);
		cvui::printf(dst, 476, 317, 0.5, 0x7fff7f, "%.4f",fftMaxX_temp);
		cvui::printf(dst, 623, 317, 0.5, 0x7fff7f, "%.4f",fftMaxY_temp);
		cvui::printf(dst, 760, 317, 0.5, 0x7fff7f, "%.4f",fftMaxMajor_temp);
		cvui::printf(dst, 896, 317, 0.5, 0x7fff7f, "%.4f",fftMaxMinor_temp);

		///////////////////////////////////更新部分///////////////////////////////////
		cvui::update();
		imshow("IGSF",dst); 
		waitKey(imgshowPeriod);
		}
	return 0;
}


void* fftWorker(void*){
	double arrFFT[12][12] = { 0 };
	vector<double> fitValue(6,0);
	//FFT画图插值部分相关定义
	cv::Mat inputMajor = cv::Mat(1,14,CV_32F);
	cv::Mat inputMinor = cv::Mat(1,14,CV_32F); 
	cv::Mat outputMajor = cv::Mat(1,SIZE/2,CV_32F);
	cv::Mat outputMinor = cv::Mat(1,SIZE/2,CV_32F); 
	//定义数组用于存储傅里叶变化的时域序列src和频域序列dst
	double fftInputX_temp[SIZE], fftInputY_temp[SIZE], fftInputMajor_temp[SIZE], fftInputMinor_temp[SIZE];
	Complex_ fftOutputX_temp[SIZE], fftOutputY_temp[SIZE], fftOutputMajor_temp[SIZE], fftOutputMinor_temp[SIZE];
	double fftOutputAmpX_temp[SIZE], fftOutputAmpY_temp[SIZE], fftOutputAmpMajor_temp[SIZE];
	double fftOutputAmpMinor_temp[SIZE], fftMaxX_temp, fftMaxY_temp, fftMaxMajor_temp, fftMaxMinor_temp;

	while(true){
		//arrFFT的初始化
		pthread_mutex_lock(&Arr.arr_lock);
		for(int i = 0; i<12; i++){
			for(int j = 0; j<12; j++){
				arrFFT[i][j] = Arr.arr[i][j];
			}
		}
		pthread_mutex_unlock(&Arr.arr_lock);
		pthread_mutex_lock(&IGSFResult.IGSFResult_lock);
		fitValue = IGSFResult.fit_value;
		pthread_mutex_unlock(&IGSFResult.IGSFResult_lock);
		///////////////////////////FFT求解过程//////////////////////////////////////////////////////////////////////////
		readSrcData(arrFFT, fftInputX_temp, fftInputY_temp, fitValue[1], fitValue[2], SIZE, isBeam);
		gauss(fftInputMajor_temp, fitValue[3]);
		gauss(fftInputMinor_temp, fitValue[4]);
		FFT_func(fftInputMajor_temp, fftOutputMajor_temp, SIZE);
		FFT_func(fftInputMinor_temp, fftOutputMinor_temp, SIZE);
		FFT_func(fftInputX_temp, fftOutputX_temp, SIZE);
		FFT_func(fftInputY_temp, fftOutputY_temp, SIZE);
		fftAmp(fftOutputX_temp, fftOutputAmpX_temp);
		fftAmp(fftOutputY_temp, fftOutputAmpY_temp);
		fftAmp(fftOutputMajor_temp, fftOutputAmpMajor_temp);
		fftAmp(fftOutputMinor_temp, fftOutputAmpMinor_temp);
		for(int i = 0; i<14 ; i++){
			inputMajor.at<float>(0,i) = fftOutputAmpMajor_temp[i];
			inputMinor.at<float>(0,i) = fftOutputAmpMinor_temp[i];
		}
		cv::resize(inputMajor, outputMajor, outputMajor.size(), 0 ,0, cv::INTER_LINEAR);
		cv::resize(inputMinor, outputMinor, outputMinor.size(), 0 ,0, cv::INTER_LINEAR);
		for(int i=0; i<SIZE/2; i++){
			fftOutputAmpMajor_temp[i] = outputMajor.at<float>(0,i);
			fftOutputAmpMinor_temp[i] = outputMinor.at<float>(0,i);
		}
		fftMaxX_temp = findMaxFreq(fftOutputAmpX_temp);
		fftMaxY_temp = findMaxFreq(fftOutputAmpY_temp);
		fftMaxMajor_temp = findMaxFreq(fftOutputAmpMajor_temp);
		fftMaxMinor_temp = findMaxFreq(fftOutputAmpMinor_temp);
		pthread_mutex_lock(&FFT.FFT_lock);
		for(int i = 0; i< SIZE; i++){
			FFT.fftOutputAmpX[i] = fftOutputAmpX_temp[i];
			FFT.fftOutputAmpY[i] = fftOutputAmpY_temp[i];
			FFT.fftOutputAmpMajor[i] = fftOutputAmpMajor_temp[i];
			FFT.fftOutputAmpMinor[i] = fftOutputAmpMinor_temp[i];
		}
		FFT.fftMaxX = fftMaxX_temp;
		FFT.fftMaxY = fftMaxY_temp;
		FFT.fftMaxMajor = fftMaxMajor_temp;
		FFT.fftMaxMinor = fftMaxMinor_temp; 
		pthread_mutex_unlock(&FFT.FFT_lock);
		usleep(fftPeriod);
	}
}

void* epicsWorker(void*){
	//初始化EPICS
	bool epicsStatus_temp = false;
	int count = 0;
	bool px3Status_temp = false;
	bool HVEnable_temp = 0;
	double hvRead_temp, integratedTimeRead_temp, px3Tempature_temp;
	double  integratedTimeSet_temp = 0;
	double px3Pressure_temp, px3Humidity_temp;
	double HVSet_temp = 0;
	double pvchannel_temp[135];

	pvID=epics_init();
	while (true)
	{
		count = 0;
		if(ca_get(DBR_ENUM, pvID[1], &px3Status_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_get(DBR_DOUBLE, pvID[2], &hvRead_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_get(DBR_DOUBLE, pvID[6], &integratedTimeRead_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_get(DBR_DOUBLE, pvID[7], &px3Tempature_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_get(DBR_DOUBLE, pvID[8], &px3Pressure_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_get(DBR_DOUBLE, pvID[9], &px3Humidity_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_array_get(DBR_DOUBLE,135,pvID[0],&pvchannel_temp) == ECA_NORMAL){
			count++;
		}
		ca_pend_io(3.0);

		pthread_mutex_lock(&PV.PV_lock);
		HVEnable_temp = PV.hvEnable;
		HVSet_temp = PV.hvSet;
		integratedTimeSet_temp = PV.integratedTimeSet;
		PV.integratedTimeRead = integratedTimeRead_temp;
		PV.px3Humidity = px3Humidity_temp;
		PV.px3Pressure = px3Pressure_temp;
		PV.px3Status = px3Status_temp;
		PV.px3Tempature = px3Tempature_temp;
		PV.hvRead = hvRead_temp;
		for(int i = 0; i<135; i++)
			PV.pvChanel[i] = pvchannel_temp[i];
		pthread_mutex_unlock(&PV.PV_lock);

		if(ca_put(DBR_ENUM, pvID[4], &HVEnable_temp) == ECA_NORMAL){
			count++;
		}
		if(HVEnable_temp && ca_put(DBR_DOUBLE, pvID[3], &HVSet_temp) == ECA_NORMAL){
			count++;
		}
		if(ca_put(DBR_DOUBLE, pvID[5], &integratedTimeSet_temp) == ECA_NORMAL){
			count++;
		}

		if(count == 10){
			myepicsStatus = true;
		}else{
			myepicsStatus = false;
		}
		usleep(epicsPeriod);
	}
}

void* igsfWorker(void*){
	const char* outFitValuePath = "../data/fit_result.txt";
	Mat image(12,12,CV_32S);
	double noiseThreshold_temp = 0.05;
	vector<double> fitValue(6,0);
	int zeroCount = 0;
	////////////////////////////IGSF求解过程////////////////////////////////////////////////////////////////////////
	//创建Preco_imgPixel类并调用成员函数对数据进行预处理
	while(true){
		//矩阵初始化
		zeroCount = 0;
		pthread_mutex_lock(&Arr.arr_lock);
		for (int i = 0; i < 12; i++){
			for (int j = 0; j < 12; j++){
				if(Arr.arr[i][j] == 0){
					zeroCount++;
				}
				image.at<int>(i, j) = int(Arr.arr[i][j]);
			}
		}
		pthread_mutex_unlock(&Arr.arr_lock);
		if(zeroCount == 144){
			usleep(10000);
			continue;
		}
		pthread_mutex_lock(&NoiseThreshold.noiseThreshold_lock);
		noiseThreshold_temp = NoiseThreshold.nosieThreshold;
		pthread_mutex_unlock(&NoiseThreshold.noiseThreshold_lock);
		Preco_imgPixel lena(12, 12, image, noiseThreshold_temp);    
		vector<int> max = lena.find_Max();                //找到像素最大值以及对应的坐标序列(i,j)
		vector<vector<double>> _B = lena.CreateCoorFit(max[2], isBeam);      //创建B矩阵
		vector<double> Fx;
		for (auto i : _B){
			Fx.push_back(i[0]);
		}
		vector<vector<double>> _A = lena.CreateA(Fx);     //创建A矩阵
		Matrix B(_B), A(_A),invS(6,6);                    //新建Matrix类 A,B并用_A,_B初始化;新建Matrix类 invS
		vector<vector<double>> _R_B = B.ReverseMat();     //对B矩阵求转置  
		Matrix R_B(_R_B);
		Matrix S = R_B.matMulti(B);                       //B的转置乘以B，得到S矩阵            
		S.inv(invS.getMat());
		Matrix TMP = invS.matMulti(R_B);                  //invR*S得到C矩阵
		Matrix C = TMP.matMulti(A);
		fitValue = C.fitParam(); 
		pthread_mutex_lock(&IGSFResult.IGSFResult_lock);
		IGSFResult.fit_value = fitValue;                         //得到6个拟合参数
		pthread_mutex_unlock(&IGSFResult.IGSFResult_lock);
		output_fitvalue(outFitValuePath, fitValue);      //输出拟合结果到文件
		usleep(igsfPeriod);
	}
}

int main(){
	pthread_t id[4];
	pthread_create(&id[0], NULL, epicsWorker, NULL);
	pthread_create(&id[1], NULL, igsfWorker, NULL);
	pthread_create(&id[2], NULL, fftWorker, NULL);
	pthread_create(&id[3], NULL, imshowWorker, NULL);
	getchar();
	return(0);
}







