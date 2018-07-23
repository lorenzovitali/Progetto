#include<iostream>
#include<vector>
#include<mpi.h>
#include<random>
#include<math.h>
#include <algorithm>

double EO_montecarlo(double S0;double sigma, double r, double T){
	
	int rank,size;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	MPI_Bcast(&S0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	MPI_Bcast(&r,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&T,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	int N = 1000;
	int n_loc = N/size;
	if(rank<N%size){
		n_loc++
	}
	
	std::default_random_engine engine (n_loc);
	std::uniform_real_distribution<double> distro (-1., 1.);
	std::vector<double> different_path(n_loc);
	
	double first_term = (r - 1/2 * sigma^2) *T ;

	
	for(int i=0;i<n_loc;i++){
		const double z = distro (engine);
 		double S = S0 * exp( first_term + sigma* sqrt(T) * z);
 		double cash_flow = max{S,0};
		different_path.pushback(cash_flow);
	}
	
	//calcolo la media locale
	double mean = 0;
	for(int i=0; i<n_loc;i++){
		mean += different_path(i);
	}
	mean = mean/n_loc;
	
	MPI_Reduce(MPI_IN_PLACE, &mean, 1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD );
	
	if(rank == 0){
		mean = mean/size;
		return mean;
	}
	
	
}
