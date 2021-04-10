#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matplotlib/matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;

const int N = 30;
const int iter = 10000;
double del_th = (10*M_PI)/180, J=1, KbT=0.01, KbTm=2,h=0;

void init_S(double** S){
	double r;	
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			r = (double)rand();
			S[i][j] = r - 2*M_PI*(int)(r/(2*M_PI));	
		}
	}
}

int mod(int i){
	if(i<0)
		return N+i;
	else
		return i%N;
}

double M_calc(double** S){
	double Mx=0, My=0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Mx += cos(S[i][j]);
			My += sin(S[j][j]);
		}
	}
	return sqrt((Mx*Mx + My*My))/(N*N);
}

double E_calc(double** S){
	double E=0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			E = -J*(cos(S[i][j]-S[mod(i+1)][j])+cos(S[i][j]-S[mod(i-1)][j])+
					cos(S[i][j]-S[i][mod(j+1)])+cos(S[i][j]-S[i][mod(j-1)]));
		}
	}

	return E;
}

void metropolis_iterate(double** S){
	double ep, en, de, r_th, p, r;
	for(int k=0; k<iter; k++){
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				ep = -h*sin(S[i][j]) - J*(cos(S[i][j]-S[mod(i+1)][j])+cos(S[i][j]-S[mod(i-1)][j])+
					cos(S[i][j]-S[i][mod(j+1)])+cos(S[i][j]-S[i][mod(j-1)]));

				r = rand();
				r_th = r - 2*del_th*(int)(r/(2*del_th)) - del_th;	

				en = -h*sin(S[i][j]+r_th) - J*(cos(S[i][j]+r_th-S[mod(i+1)][j])+cos(S[i][j]+r_th-S[mod(i-1)][j])+
					cos(S[i][j]+r_th-S[i][mod(j+1)])+cos(S[i][j]+r_th-S[i][mod(j-1)]));

				de = en-ep;
				p = exp(-de/KbT);
				r = ((double)(rand()%101))/100;
				if(r<p)
					S[i][j] += r_th;
			}
		}
	}
}

int main(){
	double** S = new double*[N];
	for(int i = 0; i < N; ++i)
    	S[i] = new double[N];

    init_S(S);
    J=1;
    h=2;
    KbT = 0.01;
    metropolis_iterate(S);

    std::vector<double> x, y, u, v;
    for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			x.push_back(i);
            u.push_back(cos(S[i][j]));
            y.push_back(j);
            v.push_back(sin(S[i][j]));
		}
	}
	plt::quiver(x, y, u, v);
    plt::show();
	return 0;
}