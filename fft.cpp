#include "fft.h"
////////////////////////////////////////////////////////////////////
//定义一个复数计算，包括乘法，加法，减法
////////////////////////////////////////////////////////////////////
void Add_Complex(Complex_ * src1,Complex_ *src2,Complex_ *dst){
    dst->imagin=src1->imagin+src2->imagin;
    dst->real=src1->real+src2->real;
}
void Sub_Complex(Complex_ * src1,Complex_ *src2,Complex_ *dst){
    dst->imagin=src1->imagin-src2->imagin;
    dst->real=src1->real-src2->real;
}
void Multy_Complex(Complex_ * src1,Complex_ *src2,Complex_ *dst){
    double r1=0.0,r2=0.0;
    double i1=0.0,i2=0.0;
    r1=src1->real;
    r2=src2->real;
    i1=src1->imagin;
    i2=src2->imagin;
    dst->imagin=r1*i2+r2*i1;
    dst->real=r1*r2-i1*i2;
}
////////////////////////////////////////////////////////////////////
//在FFT中有一个WN的n次方项，在迭代中会不断用到，具体见算法说明
///////////////////////////////////////////////////////////////////
void getWN(double n,double size_n,Complex_ * dst){
    double x=2.0*M_PI*n/size_n;
    dst->imagin=-sin(x);
    dst->real=cos(x);
}
////////////////////////////////////////////////////////////////////
//从文件中读取原始数据并对其余数据段补零
////////////////////////////////////////////////////////////////////
void readSrcData(double(&arr)[12][12], double* xData, double* yData, double xShift, double yShift, int length, bool isBeam){
    int xOffset = int(xShift/3.8);
    int yOffset = int(yShift/3.8);
    double xOddments = xShift - xOffset*3.8;
    double yOddments = yShift - yOffset*3.8;
    for(int i=0; i<length; i++){
        if(i<12){
            if(!isBeam){
                if(xOddments < 0){
                    xData[i] = arr[i][5+xOffset];
                }
                else if(xOddments == 0){
                    xData[i] = (arr[i][5+xOffset]+arr[i][6+xOffset])/2;
                }
                else{
                    xData[i] = arr[i][6+xOffset];
                }

                if(yOddments < 0){
                    yData[i] = arr[5+yOffset][i];
                }
                else if(yOddments == 0){
                    yData[i] = (arr[5+yOffset][i]+arr[6+yOffset][i])/2;   
                }
                else{
                    yData[i] = arr[6+yOffset][i];
                }
            }else{
                if(xOddments < 0){
                    xData[i] = arr[6+xOffset][i];
                }
                else if(xOddments == 0){
                    xData[i] = (arr[5+xOffset][i]+arr[6+xOffset][i])/2;
                }
                else{
                    xData[i] = arr[5+xOffset][i];
                }

                if(yOddments < 0){
                    yData[i] = arr[i][5+yOffset];
                }
                else if(yOddments == 0){
                    yData[i] = (arr[i][5+yOffset]+arr[i][6+yOffset])/2;
                }
                else{
                    yData[i] = arr[i][6+yOffset];
                }
            }
            
        }
        else{
            xData[i] = 0;
            yData[i] = 0;
        }
    }
}
////////////////////////////////////////////////////////////////////
//定义FFT的初始化数据，因为FFT的数据经过重新映射，递归结构
////////////////////////////////////////////////////////////////////
int FFT_remap(double * src,int size_n){

    if(size_n==1)
        return 0;
    double * temp=(double *)malloc(sizeof(double)*size_n);
    for(int i=0;i<size_n;i++)
        if(i%2==0)
            temp[i/2]=src[i];
        else
            temp[(size_n+i)/2]=src[i];
    for(int i=0;i<size_n;i++)
        src[i]=temp[i];
    free(temp);
    FFT_remap(src, size_n/2);
    FFT_remap(src+size_n/2, size_n/2);
    return 1;


}
////////////////////////////////////////////////////////////////////
//定义FFT，具体见算法说明，注释掉的显示部分为数据显示，可以观察结果
///////////////////////////////////////////////////////////////////
void FFT_func(double * src,Complex_ * dst,int size_n){

    FFT_remap(src, size_n);
   // for(int i=0;i<size_n;i++)
    //    printf("%lf\n",src[i]);
	//clock_t start,end;
    //start=clock();
    int k=size_n;
    int z=0;
    while (k/=2) {
        z++;
    }
    k=z;
    if(size_n!=(1<<k))
        exit(0);
    Complex_ * src_com=(Complex_*)malloc(sizeof(Complex_)*size_n);
    if(src_com==NULL)
        exit(0);
    for(int i=0;i<size_n;i++){
        src_com[i].real=src[i];
        src_com[i].imagin=0;
    }
    for(int i=0;i<k;i++){
        z=0;
        for(int j=0;j<size_n;j++){
            if((j/(1<<i))%2==1){
                Complex_ wn;
                getWN(z, size_n, &wn);
                Multy_Complex(&src_com[j], &wn,&src_com[j]);
                z+=1<<(k-i-1);
                Complex_ temp;
                int neighbour=j-(1<<(i));
                temp.real=src_com[neighbour].real;
                temp.imagin=src_com[neighbour].imagin;
                Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
                Sub_Complex(&temp, &src_com[j], &src_com[j]);
            }
            else
                z=0;
        }

    }
	for(int i=0;i<size_n;i++){
		dst[i].imagin=src_com[i].imagin;
		dst[i].real=src_com[i].real;
	}
	//end=clock();
    //printf("FFT use time :%lfs for Datasize of:%d\n",(double)(end-start)/CLOCKS_PER_SEC,size_n);

}
////////////////////////////////////////////////////////////////////
//FFT结果求幅值和归一化处理
void fftAmp(Complex_ * dst, double* amp){
    double temp, Max;
    Max = sqrt((dst[0].real)*(dst[0].real)+(dst[0].imagin)*(dst[0].imagin));
    for(int i = 0; i < SIZE; i++){
        temp = sqrt((dst[i].real)*(dst[i].real)+(dst[i].imagin)*(dst[i].imagin));
        if(temp > Max){
            Max = temp; 
        }
        amp[i] = temp;
    }
    for(int i = 0; i < SIZE; i++){
        amp[i] = amp[i]/Max;
    }
}

void gauss(double(&y)[SIZE], double sigma){
    double tmp;
    for(int i=0; i <SIZE; i++){
        y[i] = exp((-(i-SIZE/2)*(i-SIZE/2)*0.1*0.1)/(2*sigma*sigma));
    }
}

double findMaxFreq(double(&Freq)[SIZE]){
    for(int i=0; i<SIZE/2; i++ ){
        if(Freq[i]<0.01) return i*0.1316*2/SIZE;
    }
    return 0.1316;
}