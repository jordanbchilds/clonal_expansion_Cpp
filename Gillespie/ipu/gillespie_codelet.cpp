#include <cmath>
#include <cassert>
#include <poplar/Vertex.hpp>
#include <algorithm>
#include <array>
#include <ipu_vector_math>
#include <ipudef.h>
#include <ipu_builtins.h>

using namespace poplar;
using namespace std;

class sim_network_vertex : public poplar::Vertex
{
public:
	Input<int> w_init;
	Input<int> m_init;
    Input<float> reactOne_rates;
	Input<float> reactTwo_rates;
	Input<float> reactThree_rates;
	Input<float> reactFour_rates;
	Input<float> reactFive_rates;
    Input<float> conOne_rates;
	Input<float> conTwo_rates;

    Output<Vector<int>> out;
	//Output<float> out;
	
	struct sim_network {
		float Tmax;
		float step_out;
		int Nout;
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
	
	float runif_01(){
		return (float) (__builtin_ipu_urand32()) / (float) (4294967295);
	}
	
	double rand_exp(float lambda){
		float unif_01 = runif_01();
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

		float u = runif_01() ;
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
	
	float rep_controller(float con_rates[2], float rep_rate, int error) {
		float new_rate = error >= 0 ? rep_rate*2.0/(1.0+exp(error*con_rates[0])) : rep_rate*(1.0+exp(-1*error*con_rates[1]))/2.0 ;
		return new_rate;
	}

	
	void gillespied(int* x_init, float* rates, float* con_rates, int* out_array, sim_network simnet){
		int nTimes = simnet.nTimes;
		float* times = simnet.times_ptr;
		int n_species = simnet.n_species;
		int n_reactions = simnet.n_reactions;
		int* S_pt = simnet.Stoi;
		int* Pre_pt = simnet.Pre;

		int x[2];
		x[0] = *x_init; x[1] = *(x_init+1);

		int count = 0;
		float tt = 0.0;
		// float target = *times;
		int C0 = x[0]+x[1];
		int copyNum = C0;

		float temp_rates[5];
		temp_rates[2] = *(rates+2);
		temp_rates[3] = *(rates+3);
		temp_rates[4] = *(rates+4);

		while( tt <= *(times+nTimes-1) ){
			temp_rates[0] = rep_controller(con_rates, *rates, copyNum-C0);
			temp_rates[1] = rep_controller(con_rates, *(rates+1), copyNum-C0);
			float hazards[5];
			float haz_total = 0.0;
			for(int i=0; i<n_reactions; ++i){
				float h_i = temp_rates[i];
				for(int j=0; j<n_species; ++j)
					h_i *= choose(x[j], *( Pre_pt+i*n_species+j ));
				hazards[i] = h_i;
				haz_total += h_i;
			}

			tt += rand_exp(haz_total);

			while( tt >= *(times+count) && count<nTimes){
				cout<< count << " ";
				*(out_array+count*n_species) = x[0];
				*(out_array+count*n_species+1) = x[1];
				count += 1;
			}

			int r = rand_react(hazards);
			x[0]  += *( S_pt + r*n_species );
			x[1]  += *( S_pt + r*n_species + 1 );
			copyNum = x[0]+x[1];

			if(copyNum>2*C0 || copyNum==0){
				for(int i=count; i<nTimes; ++i){
					*(out_array+i*n_species) = 0;
					*(out_array+i*n_species+1) = 0;
				}
				tt = 1e99;
			}
		}
	}
	bool compute()
	{
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

		sim_network spn;
		spn.Tmax = 120.0*365.0*24.0*3600.0; // 120 years in seconds
		spn.step_out = 365.0*24.0*60.0*60.0; // 1 year in seconds
		spn.Nout = (long unsigned int) (spn.Tmax/spn.step_out+1);
		spn.n_reactions = Nreact;
		spn.n_species = Nspecies;
		spn.Post = Post_ptr;
		spn.Pre = Pre_ptr;
		spn.Stoi = S_ptr;
		
		int x_init[Nspecies];
		x_init[0] = w_init;
		x_init[1] = m_init;
		float react_rates[Nreact];
		react_rates[0] = reactOne_rates;
		react_rates[1] = reactTwo_rates;
		react_rates[2] = reactThree_rates;
		react_rates[3] = reactFour_rates;
		react_rates[4] = reactFive_rates;
		float con_rates[2];
		con_rates[0] = conOne_rates;
		con_rates[1] = conTwo_rates;
		
		int output[spn.Nout][spn.n_species];
		int* output_ptr = &output[0][0];

		gillespied(x_init, react_rates, con_rates, output_ptr, spn);

		int index = 0;
		for(int i=0; i<spn.Nout; ++i){
			out[index] = output[i][0];
			out[index+1] = output[i][1];
			index += 2;
		}
		return true;
	}
};
