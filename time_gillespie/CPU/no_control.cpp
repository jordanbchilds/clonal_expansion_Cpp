#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>

using namespace std;

struct sim_network {
	const float* times_ptr;
	float Tmax;
	int nTimes;
	int n_reactions;
	int n_species;
	int* Post;
	int* Pre;
	int* Stoi;
};

unsigned choose(unsigned n, unsigned k){
	if (k>n) return 0;
	if (k*2>n) k = n-k;
	if (k==0) return 1;
	if (k==1) return n;
	int result = n;
	for(int i=2; i<=(n-k); ++i){
		result *= (n-i+1);
		result /= i;
	}
	return result;
}

float rand_unif(float lower=0.0, float upper=1.0){
	float unif_01 = (float) rand() / (float) 4294967295 ; // 4,294,967,295 is the max value of an unsigned iteger... I think.
	return unif_01*(upper-lower) + lower;
}

double rand_exp(float lambda){
	float unif_01 = rand_unif();
	return -1.0*log(1.0-unif_01)/lambda;
}

int rand_react(float* weights){
	float norm = 0.0; // normalising constant
	float cumWeights[5];
	for(int i=0; i<5; ++i)
		norm += *(weights+i);
	
	for(int i=0; i<5; ++i){
		float cc = 0.0;
		for(int j=0; j<=i; ++j)
			cc += *(weights+j)/norm;
		cumWeights[i] = cc;
	}

	float u = rand_unif() ;
	if( 0<=u && u<cumWeights[0] ){
		return 0;
	} else {
		for(int i=1; i<5; ++i){
			if( cumWeights[i-1]<u && u<=cumWeights[i] )
				return i;
		}
	}
	return -1; // requires a return outside the loop
}

void gillespied(int* x_init, float* rates, float* con_rates, int* out_array, sim_network simnet){
	
	int nTimes = simnet.nTimes;
	float Tmax = simnet.Tmax;
	const float* times = simnet.times_ptr;
	int n_species = simnet.n_species;
	int n_reactions = simnet.n_reactions;
	int* S_pt = simnet.Stoi;
	int* Pre_pt = simnet.Pre;

	int x[2];
	x[0] = *x_init;
	x[1] = *(x_init+1);

	int count = 0;
	float tt = 0.0;
	// float target = *times;
	int C0 = x[0]+x[1];
	int copyNum = C0;

	while( tt <= Tmax ){
		float hazards[5];
		float haz_total = 0.0;
		for(int i=0; i<n_reactions; ++i){
			float h_i = *(rates+i);
			for(int j=0; j<n_species; ++j)
				h_i *= choose(x[j], *( Pre_pt+i*n_species+j ));
			hazards[i] = h_i;
			haz_total += h_i;
		}

		tt += rand_exp(haz_total);

		if( tt >= *(times+count) && count<=nTimes ){
			*(out_array+count*n_species) = x[0];
			*(out_array+count*n_species+1) = x[1];
			count += 1;
		}

		int r = rand_react(hazards);
		x[0]  += *( S_pt + r*n_species );
		x[1]  += *( S_pt + r*n_species + 1 );
		copyNum = x[0]+x[1];
		
		if(copyNum>2*C0 || copyNum==0){
			// copyNum==0 is specific to this system Darren's general equiv was haz_total<1e-10
			for(int i=count; i<nTimes; ++i){
				*(out_array+i*n_species) = 0;
				*(out_array+i*n_species+1) = 0;
			}
			tt = 1e99;
		}
	}
}

int main(){
	const int Nreact = 5;
	const int Nspecies = 2;
	int Pre_mat[Nreact][Nspecies] = { {1,0}, {0,1}, {1,0}, {0,1}, {1,0} };
	int* Pre_ptr = &Pre_mat[0][0];
	int Post_mat[Nreact][Nspecies] = { {2,0}, {0,2}, {0,0}, {0,0}, {1,1} };
	int* Post_ptr = &Post_mat[0][0];
	int S_mat[Nreact][Nspecies];
	int* S_ptr = &S_mat[0][0];
	for(int i=0; i<Nreact; ++i){
		for(int j=0; j<Nspecies; ++j)
			S_mat[i][j] = Post_mat[i][j] - Pre_mat[i][j];
	}
	
	const unsigned int Nout = 121;
	float times[Nout];
	for(int i=0; i<Nout; ++i){
		times[i] = i*365.0;
	}
	
	sim_network spn;
	spn.nTimes = Nout;
	spn.times_ptr = &times[0];
	spn.Tmax = times[Nout-1];
	spn.n_reactions = Nreact;
	spn.n_species = Nspecies;
	spn.Post = Post_ptr;
	spn.Pre = Pre_ptr;
	spn.Stoi = S_ptr;
	
	int x_init[Nspecies] = {500, 500};
	
	float react_rates[Nreact] = {2.64e-3, 2.64e-3, 2.64e-3, 2.64e-3, 0.0};
	float con_rates[2] = {2.0e-3, 2.0e-3};
	
	int output[spn.nTimes * spn.n_species];
	int* output_ptr = &output[0];
	
	const unsigned int Nsim = 1000;
	const unsigned int Ntile = 1472;
	float simTimes[Nsim] = {0};

	for(int i=0; i<Nsim; ++i){

		// chrono::time_point start = chrono::high_resolution_clock::now();
		auto start = chrono::high_resolution_clock::now();
		
		for(int j=0; j<Ntile; ++j)
			gillespied(&x_init[0], &react_rates[0], &con_rates[0], output_ptr, spn);
			
		// chrono::time_point end = chrono::high_resolution_clock::now();
		auto end = chrono::high_resolution_clock::now();
		
		chrono::duration<double, std::milli> ms_double = (end - start);
		
		std::cout << ms_double.count() << std::endl;
		simTimes[i] = ms_double.count();
	}
	__fs::filesystem::current_path("/Users/jordanchilds/Documents/GitHub/clonal_expansion_Cpp/time_gillespie/CPU");
	ofstream file ("./noControl_sim_times.txt");
	for(int j=0; j<Nsim; ++j){
		file << simTimes[j] << endl ;
	}
	file.close();
	
}
