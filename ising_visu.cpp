#include "math.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
using namespace std;
using namespace cv;

const int N = 40;
const int nd = 60;
const int iter = 1000;
double J=1, KbT=0.1, KbTm=4,h=0;

void init_S(double** S){
	double r;	
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			r = rand();			
			S[i][j] = (int)r%2+~(int)r%2;
		}
	}
}

int mod(int i){
	if(i<0)
		return N+i;
	else
		return i%N;
}

void metropolis_iterate(double** S){
	double de,p,r;
	for (int k = 0; k < iter; ++k)
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				de = 2*h*S[i][j] + 2*J*(S[i][j]*S[mod(i+1)][j]+S[i][j]*S[mod(i-1)][j]+S[i][j]*S[i][mod(j+1)]+S[i][j]*S[i][mod(j-1)]);
				p = exp(-de/KbT);
				r = ((double)((rand()+1)%101))/100;
				if(r<p)
					S[i][j] *= -1;
			}
		}
	}
}


int main(){
	// Dynamic memory allocation of NbyN matrix
	double** S = new double*[N];
	for(int i = 0; i < N; ++i)
    	S[i] = new double[N];

    init_S(S);//Initialization of the random spin configuration
   	
   	J = 0.1;
   	h = 1;
   	KbT = 1;

    metropolis_iterate(S);
    //OpenCV Mat object
    Mat M((int)N,(int)N,CV_8UC3,Vec3b(0,0,0));

    for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			//Spin up gets red BGR=(0,0,255)
			if(S[i][j]==1){
				M.at<Vec3b>(i,j)[0] = 0;
				M.at<Vec3b>(i,j)[1] = 0;
				M.at<Vec3b>(i,j)[2] = 255;
			}
			else{
				//Spin down gets yellow BGR=(0,242,250)
				M.at<Vec3b>(i,j)[0] = 0;
				M.at<Vec3b>(i,j)[1] = 242;
				M.at<Vec3b>(i,j)[2] = 250;
			}
		}
	}

	namedWindow("image",CV_WINDOW_FREERATIO);
	imshow("image",M);
	waitKey(0);
	imwrite("ising_no_field.png",M);

	return 0;
}

