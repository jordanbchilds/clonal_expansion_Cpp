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
	poplar::Input<int> w_init;
	poplar::Input<int> m_init;
    poplar::Input<float> reactOne_rates;
	poplar::Input<float> reactTwo_rates;
	poplar::Input<float> reactThree_rates;
	poplar::Input<float> reactFour_rates;
	poplar::Input<float> reactFive_rates;
    poplar::Input<float> conOne_rates;
	poplar::Input<float> conTwo_rates;
	
	poplar::InOut<poplar::Vector<int>> w_popDyn;
	poplar::InOut<poplar::Vector<int>> m_popDyn;
	
    poplar::Output<float> out;
	
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
	
	float rand_unif(float lower=0.0, float upper=1.0){
		float unif_01 = (float) (__builtin_ipu_urand32()) / (float) (4294967295) ; // 4,294,967,295 is the max value of an unsigned iteger
		float unif = unif_01*(upper-lower) + lower;
		return unif;
	}
	
	
	double rand_exp(float lambda){ // both lambda and x are positive - use type unsigned double?
		// float unif_01 = (__builtin_ipu_urand_f32()+1.0)/2.0;
		float unif_01 = rand_unif();
		return -1.0*log(1-unif_01)/lambda;
	}
	
	int rand_react(float* weights=nullptr){
		// size and weights can only be positive - use types unsigned int and unsigned double?
		float norm; // normalising constant
		if( weights!=nullptr ){ // if weights given
			norm = 0;
			for(int i=0; i<5; ++i)
				norm += *(weights+i);
		} else {
			norm = 5;
		} // if no weights given norm is the number of reactions
		
		float cumWeights[5];
		if( weights==nullptr ){
			for(int i=0; i<5; ++i)
				cumWeights[i] = (i+1)/norm;
		} else {
			for(int i=0; i<5; ++i){
				float cc = 0;
				for(int j=0; j<=i; ++j)
					cc += *(weights+j)/norm;
				cumWeights[i] = cc;
			}
		}
		float u = __builtin_ipu_urand_f32() ;
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
	}
	
	void gillespied(int* x_init, float* rates, float* con_rates, int* out_array, sim_network simnet){
		
		int n_species = simnet.n_species;
		int n_reactions = simnet.n_reactions;
		float step_out = simnet.step_out;
		float Tmax = simnet.Tmax;
		int* S = simnet.Stoi;
		int* Pre = simnet.Pre;
		
		int x[2];
		for(int i=0; i<n_species; ++i)
			x[i] = *(x_init+i);
		/*
		for(int j=0; j<n_species; ++j){
			// *(out_array+j) = x[j];
		}
		*/
		w_popDyn[0] = x[0];
		m_popDyn[0] = x[1];

		int count = 1;
		float target = step_out;
		float tt = 0;

		int C0 = x[0]+x[1];
		int copyNum = C0;
		
		while( count<simnet.Nout ){
			
			float temp_rates[5];
			temp_rates[0] = rep_controller(con_rates, *rates, copyNum-C0);
			temp_rates[1] = rep_controller(con_rates, *(rates+1), copyNum-C0);
			for(int i=2; i<n_reactions; ++i)
				temp_rates[i] = *(rates+i);
			
			float hazards[5];
			for(int i=0; i<n_reactions; ++i){
				float h_i = temp_rates[i];
				for(int j=0; j<n_species; ++j)
					h_i *= choose(x[j], *(Pre+i*n_species+j));
				hazards[i] = h_i;
			}

			float haz_total = 0;
			for(int i=0; i<n_reactions; ++i)
				haz_total += hazards[i];

			if( copyNum == 0 )
				break;
			else
			 tt += rand_exp(haz_total);
			
			if( tt>=target ){
				/*
				for(int j=0; j<n_species; ++j){
					*(out_array+count*n_species+j ) = x[j];
				}
				*/
				
				w_popDyn[count] = x[0];
				m_popDyn[count] = x[1];
				count += 1;
				target += step_out;
			}
			
			/*
			int r = rand_react(hazards);
			for(int j=0; j<n_species; ++j)
				x[j] += *(S+r*n_species+j);
			
			copyNum = x[0]+x[1];
			if( count>simnet.Nout || copyNum==0 )
				break;
			 */
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
		for(int i=0; i<Nreact; ++i){
			for(int j=0; j<Nspecies; ++j)
				S_mat[i][j] = Post_mat[i][j] - Pre_mat[i][j];
		}
		int* S_ptr = &S_mat[0][0];

		sim_network spn;
		spn.Tmax = 365.0*24.0*3600.0; // 1 year in seconds
		spn.step_out = 10.0*24.0*60.0*60.0; // ten days in seconds
		spn.Nout = (long unsigned int) (spn.Tmax/spn.step_out + 1.0);
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

		*out = w_popDyn[10]+m_popDyn[10];
		return true;
	}
};
