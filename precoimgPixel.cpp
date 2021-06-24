#include"precoimgPixel.h"



Preco_imgPixel::Preco_imgPixel(int row, int col, Mat& image, double w):imagePixel(image) {
	weight = w;
}


void Preco_imgPixel::showImage() {
	if (imagePixel.empty()) // Check for invalid input
	{
		cout << "Could not open or find the image" << std::endl;
		return;
	}


	namedWindow("Display window", WINDOW_AUTOSIZE); // Create a window for display.
	imshow("Display window", imagePixel); // Show our image inside it.

	waitKey(0); // Wait for a keystroke in the window
}

vector<int> Preco_imgPixel::find_Max() {
	vector<int> max(3,0);
	max[0] = 0;
	max[1] = 0;
	max[2]= int(imagePixel.at<int>(0, 0));
	for (int i = 0; i < imagePixel.rows; i++)
		for (int j = 0; j < imagePixel.cols; j++)
		{
			if (max[2] < int(imagePixel.at<int>(i, j)))
			{
				max[0] = i;
				max[1] = j;
				max[2] = imagePixel.at<int>(i, j);
			}
		}
	return max;
}



vector<vector<double>> Preco_imgPixel::CreateCoorFit(short max, bool isBeam) {
	vector<vector<double>> Coor_fit;
	vector<double> tmp(3,0);
	vector<double> rows(6, 0);
	for (int i = 0; i < imagePixel.rows; i++)
		for (int j = 0; j < imagePixel.cols; j++)
		{
			if (int(imagePixel.at<int>(i, j))>=max*weight)
			{
				if(isBeam){
					tmp[0] = j* PIX_SIZE - PIX_SIZE*5.5;
					tmp[1] = (11-i) * PIX_SIZE - PIX_SIZE * 5.5;
					tmp[2] = int(imagePixel.at<int>(i, j));
				}
				else{
					tmp[0] = i* PIX_SIZE - PIX_SIZE*5.5;
					tmp[1] = j * PIX_SIZE - PIX_SIZE * 5.5;
					tmp[2] = int(imagePixel.at<int>(i, j));
				}

				rows[0] = tmp[2];
				rows[1] = tmp[0] * tmp[2];
				rows[2] = tmp[1] * tmp[2];
				rows[3] = tmp[0] * tmp[1] * tmp[2];
				rows[4] = tmp[2] * tmp[0] * tmp[0];
				rows[5] = tmp[2] * tmp[1] * tmp[1];

				Coor_fit.push_back(rows);
			}
		}
	return Coor_fit;
}


vector<vector<double>> Preco_imgPixel::CreateA(vector<double>& Fx) {
	vector<vector<double>> res;
	vector<double> tmp(1,0);
	for (int i = 0; i < Fx.size(); i++) {
		tmp[0] = Fx[i] * log(Fx[i]);
		res.push_back(tmp);
	}
	return res;
}


void coor_Convert(const char* inpath, double(&arr)[12][12],double theta) {
	double t1, t2;
	int x = 0;
	int y = 0;
	ifstream in;
	in.open(inpath, ios::in);   //打开基于仿真数据的粒子坐标文件（注意修改文件路径）
	if (!in.is_open()) {
		cout << "Open proton.txt file failure" << endl;
	}
	int i = 0;//文件打开出错处理
	while (i++< PARTIClE_NUM) {
		in >> t1 >> t2;
		// double t1_1 = t1* cos(theta) -t2* sin(theta);
		// double t2_1 = t1 * sin(theta) + t2* cos(theta);          //旋转坐标变换
		double t1_1 = t1 * cos(theta) + t2* sin(theta);
		double t2_1 = t2* cos(theta) - t1* sin(theta);
		t1_1 += PIX_SIZE * 6;
		t2_1 += PIX_SIZE * 6;                                            //坐标变换，将粒子坐标向左上方移动（3.8*6），最下方的像素点(3.8*3.8)对应于arr[0][0]
		if (t1_1 >= 0 && t1_1 <= PIX_SIZE * 12 && t2_1 >= 0 && t2_1 <= PIX_SIZE * 12) {  //判断粒子坐标是否在像素电离室探测范围内   
			x = (int)(t1_1 / PIX_SIZE);
			y = (int)(t2_1 / PIX_SIZE);   //将粒子坐标变换到对应的像素点编号(x,y)
			arr[x][y]++;         //统计每个像素点编号表示范围内的粒子数量
		}
	}
	in.close();                     //关闭文件   
}


void output_to_File(const char* outpath, double(&arr)[12][12]) {
	ofstream out;
	out.open(outpath, ios::app);
	if (!out.is_open())
		cout << "Create file failure!" << endl;
	out<<"new data--------------------------------------------"<<endl;
	for (int i = 0; i < 12; i++) {        //将统计结果写入到目标文件中
		for (int j = 0; j < 12; j++) {
			out << arr[i][j] << " ";
		}
		out << endl;
	}
	out.close();
}

vector<chid> epics_init(){
	const char* pvName1 = "r_i128_in_current";
	const char* pvNmae2 = "r_i128_in_connection_status";
	const char* pvName3 = "r_i128_in_external_hv";
	const char* pvName4 = "c_i128_out_hv";
	const char* pvName5 = "c_i128_out_enable_external_hv";
	const char* pvName6 = "c_i128_out_integration_time";
	const char* pvName7 = "r_i128_in_integration_time";
	const char* pvName8 = "r_i128_in_ic_temp";
	const char* pvName9 = "r_i128_in_ic_pressure";
	const char* pvName10 = "r_i128_in_ic_humidty";
	vector<chid> pvID(10);
	ca_task_initialize();

	ca_search(pvName1, &pvID[0]);
	ca_search(pvNmae2, &pvID[1]);
	ca_search(pvName3, &pvID[2]);
	ca_search(pvName4, &pvID[3]);
	ca_search(pvName5, &pvID[4]);
	ca_search(pvName6, &pvID[5]);
	ca_search(pvName7, &pvID[6]);
	ca_search(pvName8, &pvID[7]);
	ca_search(pvName9, &pvID[8]);
	ca_search(pvName10, &pvID[9]);
	ca_pend_io(5.0);
	return pvID;
}


int channel_to_Pixel(chid pvID, double (&pChanel)[135], double(&arr)[12][12]){

	int status;
	status = ca_array_get(DBR_DOUBLE,135,pvID,&pChanel);
	ca_pend_io(3.0);

	for(int i=0; i<12; i++)
		for(int j=0; j<12; j++){
			arr[i][j]=0;
		}

	arr[5][6]=pChanel[0];
	arr[5][7]=pChanel[1];
	arr[5][8]=pChanel[2];
	arr[5][9]=pChanel[3];
	arr[5][10]=pChanel[4];
	arr[5][11]=pChanel[5];
	arr[4][6]=pChanel[6];
	arr[4][7]=pChanel[7];
	arr[4][8]=pChanel[8];
	arr[4][9]=pChanel[9];
	arr[4][10]=pChanel[10];
	arr[4][11]=pChanel[11];
	arr[3][6]=pChanel[12];
	arr[3][7]=pChanel[13];
	arr[3][8]=pChanel[14];
	arr[3][9]=pChanel[15];
	arr[3][10]=pChanel[16];
	arr[3][11]=pChanel[17];
	arr[2][6]=pChanel[18];
	arr[2][7]=pChanel[19];
	arr[2][8]=pChanel[20];
	arr[2][9]=pChanel[21];
	arr[2][10]=pChanel[22];
	arr[1][6]=pChanel[23];
	arr[1][7]=pChanel[24];
	arr[1][8]=pChanel[25];
	arr[1][9]=pChanel[26];
	arr[0][6]=pChanel[27];
	arr[0][7]=pChanel[28];
	arr[0][8]=pChanel[29];
	arr[6][6]=pChanel[96];
	arr[7][6]=pChanel[97];
	arr[8][6]=pChanel[98];
	arr[9][6]=pChanel[99];
	arr[10][6]=pChanel[100];
	arr[11][6]=pChanel[101];
	arr[6][7]=pChanel[102];
	arr[7][7]=pChanel[103];
	arr[8][7]=pChanel[104];
	arr[9][7]=pChanel[105];
	arr[10][7]=pChanel[106];
	arr[11][7]=pChanel[107];
	arr[6][8]=pChanel[108];
	arr[7][8]=pChanel[109];
	arr[8][8]=pChanel[110];
	arr[9][8]=pChanel[111];
	arr[10][8]=pChanel[112];
	arr[11][8]=pChanel[113];
	arr[6][9]=pChanel[114];
	arr[7][9]=pChanel[115];
	arr[8][9]=pChanel[116];
	arr[9][9]=pChanel[117];
	arr[10][9]=pChanel[118];
	arr[6][10]=pChanel[119];
	arr[7][10]=pChanel[120];
	arr[8][10]=pChanel[121];
	arr[9][10]=pChanel[122];
	arr[6][11]=pChanel[123];
	arr[7][11]=pChanel[124];
	arr[8][11]=pChanel[125];
	arr[6][5]=pChanel[64];
	arr[6][4]=pChanel[65];
	arr[6][3]=pChanel[66];
	arr[6][2]=pChanel[67];
	arr[6][1]=pChanel[68];
	arr[6][0]=pChanel[69];
	arr[7][5]=pChanel[70];
	arr[7][4]=pChanel[71];
	arr[7][3]=pChanel[72];
	arr[7][2]=pChanel[73];
	arr[7][1]=pChanel[74];
	arr[7][0]=pChanel[75];
	arr[8][5]=pChanel[76];
	arr[8][4]=pChanel[77];
	arr[8][3]=pChanel[78];
	arr[8][2]=pChanel[79];
	arr[8][1]=pChanel[80];
	arr[8][0]=pChanel[81];
	arr[9][5]=pChanel[82];
	arr[9][4]=pChanel[83];
	arr[9][3]=pChanel[84];
	arr[9][2]=pChanel[85];
	arr[9][1]=pChanel[86];
	arr[10][5]=pChanel[87];
	arr[10][4]=pChanel[88];
	arr[10][3]=pChanel[89];
	arr[10][2]=pChanel[90];
	arr[11][5]=pChanel[91];
	arr[11][4]=pChanel[92];
	arr[11][3]=pChanel[93];
	arr[5][5]=pChanel[32];
	arr[4][5]=pChanel[33];
	arr[3][5]=pChanel[34];
	arr[2][5]=pChanel[35];
	arr[1][5]=pChanel[36];
	arr[0][5]=pChanel[37];
	arr[5][4]=pChanel[38];
	arr[4][4]=pChanel[39];
	arr[3][4]=pChanel[40];
	arr[2][4]=pChanel[41];
	arr[1][4]=pChanel[42];
	arr[0][4]=pChanel[43];
	arr[5][3]=pChanel[44];
	arr[4][3]=pChanel[45];
	arr[3][3]=pChanel[46];
	arr[2][3]=pChanel[47];
	arr[1][3]=pChanel[48];
	arr[0][3]=pChanel[49];
	arr[5][2]=pChanel[50];
	arr[4][2]=pChanel[51];
	arr[3][2]=pChanel[52];
	arr[2][2]=pChanel[53];
	arr[1][2]=pChanel[54];
	arr[5][1]=pChanel[55];
	arr[4][1]=pChanel[56];
	arr[3][1]=pChanel[57];
	arr[2][1]=pChanel[58];
	arr[5][0]=pChanel[59];
	arr[4][0]=pChanel[60];
	arr[3][0]=pChanel[61];
	return status;
}


double find_arr_max_double(double (&arr)[12][12]){
	double max = arr[0][0];
	for(int i=0; i<12; i++)
		for(int j=0; j<12; j++){
			if(arr[i][j] > max){
				max = arr[i][j];
			}
		}
	return max;
}

void Pixel_to_GrayScale(double maxChanel, double (&arr)[12][12]){
	for(int i=0; i<12; i++)
		for(int j=0; j<12; j++){
			if(arr[i][j]<0){
				arr[i][j] = -arr[i][j];
			}
			else{
				arr[i][j] = arr[i][j]/maxChanel*254;
			}
		}
}

void output_fitvalue(const char* outpath, vector<double> &fitValue){
	ofstream out;
	struct tm* timenow;
	out.open(outpath, ios::app);
	if(!out.is_open()){
		cout<<"open outFitValuePath failed"<<endl;
	}
	time_t second = time(0);
	timenow = localtime(&second);
	out<<timenow->tm_year+1900<<"."<<timenow->tm_mon+1<<"."<<timenow->tm_mday<<"--"<<timenow->tm_hour<<":"<<timenow->tm_min<<":"<<timenow->tm_sec<<"		";
	for (auto i : fitValue)
		out << i << "       ";
	out << endl;
}

double AverageRandom(double min, double max)//生成一个平均分布的随机数
{
	int minInteger = (int)(min * 10000);
	int maxInteger = (int)(max * 10000);
	int randInteger = rand()*rand();
	int diffInteger = maxInteger - minInteger;
	int resultInteger = randInteger % diffInteger + minInteger;
	return resultInteger / 10000.0;
}
double Normal2(double x, double y, double miu1, double miu2, Mat matrix, float detCov) //二维正态分布密度函数
{
	float a = matrix.at<float>(0, 0);
	float b = matrix.at<float>(0, 1);
	float c = matrix.at<float>(1, 0);
	float d = matrix.at<float>(1, 1);
 
	double A = 2 * 3.1415926*detCov;
	double B = a*(x - miu1)*(x - miu1) + (b + c)*(x - miu1)*(y - miu2) + d*(y - miu2)*(y - miu2);
 
	//cout << "==" << 1.0 / A*exp(-0.5*B)<<endl;
	return 1.0/A*exp(-0.5*B);
 
}
 
Point2i NormalRandom2(double miu1, double miu2, Mat matrix, float detCov, double min, double max)//产生二维正态分布随机数
{
	Point2i p;
	double z;
	double dScope;
	do
	{
		p.x = (int)AverageRandom(min, max);
		p.y = (int)AverageRandom(min, max);
		z = Normal2(p.x, p.y, miu1, miu2, matrix, detCov);
		dScope = AverageRandom(0, Normal2(miu1, miu2, miu1, miu2, matrix, detCov));
		//cout << "*******" << dScope << "******" << endl;
	} while (dScope >= z);
	return p;
}//产生二维正态分布随机数
 
void covMatrix(float a, float b, float c, float d, Mat& inverse_matrix, float& detCov)  //协方差矩阵
{
	Mat cov_matrix(2, 2, CV_32FC1);    //协方差矩阵
	cov_matrix.at<float>(0, 0) = a;
	cov_matrix.at<float>(0, 1) = b;
	cov_matrix.at<float>(1, 0) = c;
	cov_matrix.at<float>(1, 1) = d;
	detCov = abs(a*d - b*c);
 
	//Mat inverse_matrix(2, 2, CV_32FC1);   //协方差逆矩阵
	inverse_matrix.at<float>(0, 0) = cov_matrix.at<float>(1, 1) / detCov;
	inverse_matrix.at<float>(0, 1) = -cov_matrix.at<float>(0, 1) / detCov;
	inverse_matrix.at<float>(1, 0) = -cov_matrix.at<float>(1, 0) / detCov;
	inverse_matrix.at<float>(1, 1) = cov_matrix.at<float>(0, 0) / detCov;
 
}


void coor_Convert_Sim_Version(double(&arr)[12][12],double theta) {
	double t1, t2;
	int x = 0;
	int y = 0;
	int i = 0;
	float detCov = 36;   //协方差矩阵det
	Mat inverse_matrix(2, 2, CV_32FC1);   //协方差逆矩阵
	covMatrix(6, 0, 0, 6, inverse_matrix, detCov);
	int min = 0;
	int max = 200;
	Point2i p;
	while (i++< PARTIClE_NUM) {
		p = NormalRandom2(0, 0, inverse_matrix,detCov,min,max);
		t1 = p.x;
		t2 = p.y;
		// double t1_1 = t1* cos(theta) -t2* sin(theta);
		// double t2_1 = t1 * sin(theta) + t2* cos(theta);          //旋转坐标变换
		double t1_1 = t1 * cos(theta) + t2* sin(theta);
		double t2_1 = t2* cos(theta) - t1* sin(theta);
		t1_1 += PIX_SIZE * 6;
		t2_1 += PIX_SIZE * 6;                                            //坐标变换，将粒子坐标向左上方移动（3.8*6），最下方的像素点(3.8*3.8)对应于arr[0][0]
		if (t1_1 >= 0 && t1_1 <= PIX_SIZE * 12 && t2_1 >= 0 && t2_1 <= PIX_SIZE * 12) {  //判断粒子坐标是否在像素电离室探测范围内   
			x = (int)(t1_1 / PIX_SIZE);
			y = (int)(t2_1 / PIX_SIZE);   //将粒子坐标变换到对应的像素点编号(x,y)
			arr[x][y]++;         //统计每个像素点编号表示范围内的粒子数量
		}
	}
}
