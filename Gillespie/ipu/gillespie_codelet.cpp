#include <cmath>
#include <cassert>
#include <poplar/Vertex.hpp>
#include <algorithm>
#include <array>
#include <ipu_vector_math>

using namespace poplar;
using namespace std;

class sim_network_vertex : public poplar::Vertex{
	poplar::Input<int> x_init[2];
    poplar::Input<float> react_rates[5];
    poplar::Input<float> con_rates[2];

    poplar::Output<float> out;

bool compute(){
	
class sim_network {
public:
	float Tmax;
	float step_out;
	int Nout;
	int n_reactions;
	int n_species;
	int* Post;
	int* Pre;
	int* Stoi;

	void set_mats(int* post_ptr, int* pre_ptr, int* stoi_ptr );
	sim_network(int nReacts, int nSpecies, int* post_ptr, int* pre_ptr, int* stoi_ptr, float tmax, float stepOut);
	sim_network();
};
sim_network::sim_network(){
	Tmax = 0.0;
	step_out = 0.0;
	n_reactions = 0;
	n_species = 0;
	Post = nullptr;
	Pre = nullptr;
	Stoi = nullptr;
}

sim_network::sim_network(int nReacts, int nSpecies, int* post_ptr, int* pre_ptr, int* stoi_ptr, float tmax, float stepOut){
	Tmax = tmax;
	step_out = stepOut;
	n_reactions = nReacts;
	n_species = nSpecies;
	Post = post_ptr;
	Pre = pre_ptr;
	Stoi = stoi_ptr;
	
	for(int i=0; i<n_reactions; ++i){
		for(int j=0; j<n_species; ++j){
			*(Stoi+i*nSpecies+j) = *(Post+i*nSpecies+j) - *(Pre+i*nSpecies+j);
		}
	}
}
void sim_network::set_mats(int* post_ptr, int* pre_ptr, int* stoi_ptr ){
	assert(post_ptr!=nullptr & post_ptr!=nullptr & n_reactions!=0 & n_species!=0);
	Post = post_ptr;
	Pre = pre_ptr;
	Stoi = stoi_ptr;
	
	for(int i=0; i<n_reactions; ++i){
		for(int j=0; j<n_species; ++j){
			*(Stoi+i*n_species+j) = *(Post+i*n_species+j) - *(Pre+i*n_species+j);
		}
	}
}

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

double rand_unif(float lower=0.0, float upper=1.0){
	float unif_01;
    float unif;
    unif_01 = (float) rand() / (float) RAND_MAX;
    unif = unif_01*(upper-lower) + lower;
    return unif;
}

double rand_exp(float lambda){ // both lambda and x are positive - use type unsigned double?
    float x;
    float unif_01;
    unif_01 = drand48();
    x = -log(1-unif_01)/lambda;
    return x;
}

int rand_disc(int size, float* weights=0){
    // size and weights can only be positive - use types unsigned int and unsigned double?
    float norm; // normalising constant
    if( weights!=0 ){ // if weights given
        norm = 0;
        for(int i=0; i<size; ++i)
            norm += weights[i];
    } else {
        norm = size;
    } // if no weights given norm is the
    
    float cumWeights[size];
    if( weights==0 ){
        for(int i=0; i<size; ++i)
            cumWeights[i] = (i+1)/norm;
    } else {
        for(int i=0; i<size; ++i){
            float cc = 0;
            for(int j=0; j<=i; ++j)
                cc += weights[j]/norm;
            cumWeights[i] = cc;
        }
    }
    float u = drand48();
    if( u<cumWeights[0] ){
        return 0;
    } else {
        for(int i=1; i<size; ++i){
            if( cumWeights[i-1]<u & u<=cumWeights[i] )
				return i;
        }
    }
    return -1; // requires a return outside the loop
}

float rep_controller(float con_rates[2], float rep_rate, int error) {
	if( error >= 0){
		return rep_rate*2.0/(1.0+exp(error*con_rates[0]));
	} else {
		return rep_rate*(1.0+exp(-1*error*con_rates[1]))/2.0;
	}
}

void myHazard(float* haz_ptr, int* x, float* con_rates, float* rates, int error, sim_network simnet){
	int n_reactions = simnet.n_reactions;
	int n_species = simnet.n_species;
	int* Pre = simnet.Pre;
	
	*haz_ptr = rep_controller(con_rates, *rates, error);
	*(haz_ptr+1) = rep_controller(con_rates, *(rates+1), error);
	for(int i=2; i<n_reactions; ++i)
		*(haz_ptr+i) = *(rates+i);

	for(int i=0; i<n_reactions; ++i){
		for(int j=0; j<n_species; ++j)
			*(haz_ptr+i) *= choose(x[j], *(Pre+i*n_species+j));
	}
	// return hazards;
}

void gillespied(int* x_init, float* rates, float* con_rates, int* out_array, sim_network simnet){
	// could use pointers and such instead of stating array size in function argument?
	int x[2];
	for(int i=0; i<2; ++i)
		x[i] = *(x_init+i);
	
	int n_species = simnet.n_species;
	int n_reactions = simnet.n_reactions;
	float step_out = simnet.step_out;
	float Tmax = simnet.Tmax;
	int* S = simnet.Stoi;
	int* Pre = simnet.Pre;
	
    for(int j=0; j<n_species; ++j)
		*(out_array+j) = x[j];
    
    int count = 1;
    float target = step_out;
    float tt = 0;

	int C0 = x[0]+x[1];
	int copyNum = C0;
	
    while( tt<Tmax ){
		float temp_rates[n_reactions];
		temp_rates[0] = rep_controller(con_rates, *rates, copyNum-C0);
		temp_rates[1] = rep_controller(con_rates, *(rates+1), copyNum-C0);
		for(int i=2; i<n_reactions; ++i)
			temp_rates[i] = *(rates+i);
		
		float hazards[n_reactions];
		for(int i=0; i<n_reactions; ++i){
			float h_i = temp_rates[i];
			for(int j=0; j<n_species; ++j)
				h_i *= choose(x[j], *(Pre+i*n_species+j));
			hazards[i] = h_i;
		}
		
		//float hazards[n_reactions];
		//myHazard(hazards, x, con_rates, rates, C0-copyNum, simnet);
		
        float haz_total = 0;
        for(int i=0; i<n_reactions; ++i)
			haz_total += hazards[i];
        
        if( copyNum == 0 )
			break;
        else
            tt += rand_exp(haz_total);
		if( tt>=target ){
			for(int j=0; j<n_species; ++j)
				*(out_array+count*n_species+j ) = x[j];
			count += 1;
			target += step_out;
		}
		
        int r = rand_disc(n_reactions, hazards);
		
        for(int j=0; j<n_species; ++j)
			x[j] += *(S+r*n_species+j);
		
		copyNum = x[0]+x[1];
		if( count>(simnet.Tmax/simnet.step_out + 1.0) | copyNum==0 )
			break;
    }
    // return out_arrAY;
}
	
	const int Nreact = 5;
	const int Nspecies = 2;
	int Pre_mat[Nreact][Nspecies] = { {1,0}, {0,1}, {1,0}, {0,1}, {1,0} };
	int* Pre_ptr = &Pre_mat[0][0];
	int Post_mat[Nreact][Nspecies] = { {2,0}, {0,2}, {0,0}, {0,0}, {1,1} };
	int* Post_ptr = &Post_mat[0][0];
	int S_mat[Nreact][Nspecies];
	int* S_ptr = &S_mat[0][0];
	
    //int x_init[2] = {500,500};
    float tmax = 3784320000.0; //120*365*24*60*60 in seconds
    float stepOut = 365.0*24.0*60.0*60.0; // in seconds
    //float react_rates[5] = { 3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0.0};
	//float con_rates[2] = {2.0e-3, 2.0e-3};

	
	sim_network spn = sim_network(Nreact, Nspecies, Post_ptr, Pre_ptr, S_ptr, tmax, stepOut);
	
	int output[int(spn.Tmax/spn.step_out + 1.0)][spn.n_species];
	int* output_ptr = &output[0][0];
	
	srand((unsigned)time(NULL));
	gillespied(x_init, react_rates, con_rates, output_ptr, spn);
	*out = output;
    return true;
}
}