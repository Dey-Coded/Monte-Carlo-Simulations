#include "math.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "mpi.h"
#include "matplotlib/matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;

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
			S[i][j] = (int)r%2+~(int)r%2; // 1 or -1 is returned randomly
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
				r = ((double)((rand()+1)%101))/100; // a random number in range [0,1] is returned
				if(r<p)
					S[i][j] *= -1; //spin is thus flipped with a probability of p
			}
		}
	}
}

double M_calc(double** S){
	double M=0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			M += S[i][j];
		}
	}
	return fabs(M)/double(N*N);
}

double E_calc(double** S){
	double E=0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			E -= h*S[i][j] + J*(S[i][j]*S[mod(i+1)][j]+S[i][j]*S[mod(i-1)][j]+S[i][j]*S[i][mod(j+1)]+S[i][j]*S[i][mod(j-1)]);
		}
	}
	return E;
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

	init_S(S);	//Initialization of the random spin configuration
	double dKbT = (KbTm-KbT)/((double)(nd-1)); //Seperation between KbT(hence temperature) values at which simulations will be done
	std::vector<double> x,M,E; //The vectors to store KbT, average magentization and total energy
	
	double* data = new double[3]; // (KbT,<M>,E) data that will be computed in every core(process)
	double* data_r = new double[3];// data that will be recieved by process(core) 0 by 1,2 and 3
	data[0] = KbT + dKbT*((double)rank);// initial value of KbT for every process

	do{
		i++;
		KbT = data[0];
		metropolis_iterate(S); //Running the Metropolis iterations in all processes starting with the same initial S
		//Note that in the magnetization plot generated, at low temperatures sets of 4 dots are parallel
		//This is because in the region of order(low KbT), if all the processes start their iteration with the same initial S
		//then for slightly different KbTs the final configurations are approximately the same.
		//This is not the case for temperatures above the Curie point.
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
	}while(i<nd/np);//If 60 simulations are to be done in 4 cores, the loop runs for 60/4 = 15 times

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
