#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matplotlib/matplotlibcpp.h"
#include "mpi.h"
namespace plt = matplotlibcpp;
using namespace std;

const int N = 30;
const int nd = 40;
const int iter = 10000;
double del_th = (10*M_PI)/180, J=1, KbT=0.01, KbTm=2,h=0;

void init_S(double** S){
	double r;	
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			r = (double)rand();
			S[i][j] = r - 2*M_PI*(int)(r/(2*M_PI));	// 1 or -1 is returned randomly
		}
	}
}
//Use of mod() imposes periodic boundary conditions
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
				ep = - J*(cos(S[i][j]-S[mod(i+1)][j])+cos(S[i][j]-S[mod(i-1)][j])+
					cos(S[i][j]-S[i][mod(j+1)])+cos(S[i][j]-S[i][mod(j-1)]));

				r = rand();
				r_th = r - 2*del_th*(int)(r/(2*del_th)) - del_th;	// a random number in range [-del_th,del_th] is returned

				en = - J*(cos(S[i][j]+r_th-S[mod(i+1)][j])+cos(S[i][j]+r_th-S[mod(i-1)][j])+
					cos(S[i][j]+r_th-S[i][mod(j+1)])+cos(S[i][j]+r_th-S[i][mod(j-1)]));

				de = en-ep;
				p = exp(-de/KbT);
				r = ((double)(rand()%101))/100;
				if(r<p)
					S[i][j] += r_th; //spin is thus rotated by r_th with a probability of p
			}
		}
	}
}

int main(){
	int i=0,np,rank,err;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Dynamic memory allocation of NbyN matrix
	double** S = new double*[N];
	for(int i = 0; i < N; ++i)
    	S[i] = new double[N];

	init_S(S); //Initialization of the random spin configuration
	double dKbT = (KbTm-KbT)/((double)(nd-1)); //Seperation between KbT(hence temperature) values at which simulations will be done
	std::vector<double> x,M,E; //The vectors to store KbT, average magentization and total energy
	
	double* data = new double[3]; // (KbT,<M>,E) data that will be computed in every core(process)
	double* data_r = new double[3];// data that will be recieved by process(core) 0 by 1,2 and 3
	data[0] = KbT + dKbT*((double)rank);// initial value of KbT for every process
	do{
		i++;
		KbT = data[0];
		metropolis_iterate(S);//Running the Metropolis iterations in all processes starting with the same initial S
		data[1] = M_calc(S);
		data[2] = E_calc(S);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank>0){
			MPI_Send(data,3,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		}
		else{
			for (int j = 1; j < 4; ++j)
			{
				MPI_Status status;
				MPI_Probe(j, j, MPI_COMM_WORLD, &status);
				MPI_Recv(data_r,3,MPI_DOUBLE,j,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				x.push_back(data_r[0]);
				M.push_back(data_r[1]);
				E.push_back(data_r[2]);
				cout<<"KbT: "<<data_r[0]<<" rank: "<<j<<endl;
			}
			x.push_back(data[0]);
			M.push_back(data[1]);
			E.push_back(data[2]);
			cout<<"KbT: "<<data[0]<<" rank: 0"<<endl<<endl;	
			}
		data[0] += dKbT*((double)np);

		MPI_Barrier(MPI_COMM_WORLD);
	}while(i<nd/np); //If 40 simulations are to be done in 4 cores, the loop runs for 40/4 = 10 times

	if(rank==0){
		plt::title("Equilibrium magnetization");
		plt::ylabel("<M>/N^2");
		plt::xlabel("KbT/J");
		plt::scatter(x,M);
		plt::xlim(0.0,KbTm);
		plt::show();

		plt::title("Equilibrium energy");
		plt::ylabel("<E>");
		plt::xlabel("KbT/J");
		plt::scatter(x,E);
		plt::xlim(0.0,KbTm);
		plt::show();
	}
	MPI_Finalize();

	return 0;
}