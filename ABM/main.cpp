#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace std;

// global variables for mtDNA model
const int MAX_LEN = 2000;

float rand_unif(float lower=0, float upper=1){
	srand( (unsigned)time(NULL) );
	float unif_01;
	float unif;
	unif_01 = (float) rand() / (float) RAND_MAX;
	unif = unif_01*(upper-lower) + lower;
	return unif;
}


class agent {
private:
	int unique_id;
	int parent_id;
	int species; // wild type: 0, mutant type: 1
	int n_reactions;
	float* rates;
	
public:
	agent(int id, int pid, int mol_type, float* react_rates);
	agent();
	~agent();
	
	void set_id(int id) {
		unique_id = id;
	}
	int get_id() {
		return unique_id;
	}
	void set_pid(int pid) {
		parent_id = pid;
	}
	int get_pid() {
		return parent_id;
	}
	void set_type(int mol_type) {
		assert( mol_type==0 || mol_type==1);
		species = mol_type;
		if(species==0){
			n_reactions = 3;
		} else {
			n_reactions = 3;
		}
	}
	int get_type() {
		return species;
	}
	float* get_rates() {
		return rates;
	}
	int get_nreact() {
		return n_reactions;
	}
	int react(float deltaT, float con_rates[2], float error=0);
};

// constructor
agent::agent(){
	unique_id = -1;
	parent_id = -1;
	species = -1;
	n_reactions = -1;
	rates = nullptr;
}

agent::agent(int id, int pid, int mol_type, float* react_rates){
	assert(mol_type==0 || mol_type==1);
	unique_id = id;
	parent_id = pid;
	species = mol_type;
	rates = react_rates;
	if( mol_type==0 ){
		n_reactions = 3;
	} else {
		n_reactions = 3;
	}
}
// deconstructor
agent::~agent(){}

float rep_controller(float con_rates[2], float rep_rate, float error) {
	if( error>=0){
		float out = rep_rate*2.0 / (1.0+exp(error*con_rates[0]));
		return out;
	} else {
		float out = rep_rate*(1.0+exp(-1*error*con_rates[1]))/2.0;
		return out;
	}
}

/*
float rep_controller(float con_rates[2], float rep_rate, float error){
	return 3.06e-8;
}
*/
int agent::react(float deltaT, float con_rates[2], float error) {
	/*
	 returns a integer [0,1,2,3] depending on which reaction takes place
	 0: no reaction
	 1: replication
	 2: degradation
	 3: mutation (can only occur for wild type)
	 */
	if( species==-1 )
		return -1;
	
	float temp_rates[n_reactions];
	temp_rates[0] = rep_controller(con_rates, *rates, error);
	for(int i=1; i<n_reactions; ++i)
		temp_rates[i] = *(rates+i);
	
	float haz_total = 0;
	for(int i=0; i<n_reactions; ++i)
		haz_total += temp_rates[i];
	if( haz_total==0 )
		return 0;
	

	float event_prob = 1 - exp(-1*haz_total*deltaT);
	float u1 = drand48();

	if( u1<event_prob ){
		// could change this to use the rand_disc function you wrote? - in Gillespie file
		float cumProb[n_reactions];
		cumProb[0] = temp_rates[0] /haz_total;
		
		for(int i=1; i<n_reactions; ++i){
			float p_i = temp_rates[0] /haz_total;
			for(int j=1; j<=i; ++j)
				p_i += temp_rates[j] /haz_total;
			cumProb[i] = p_i;
		}
		float u2 = drand48();
		
		if( u2<=cumProb[0] )
			return 1; // replication: 1
		for(int i=1; i<n_reactions; ++i){
			if(cumProb[i-1]<u2 & u2<=cumProb[i])
				return i+1; // degradation: 2, mutation: 3
		}
	}
	return 0; // no reaction: 0
}

void system_update(agent* system_ptr, int copyNum, int error, float step, int current_id, float* wld_rates, float* mnt_rates, float* con_rates){
	
	int react_vector[copyNum];
	
	// pararallellelise
	agent old_system[copyNum];
	for( int i=0; i<copyNum; ++i){
		old_system[i] = *(system_ptr+i);
		react_vector[i] = ( *(system_ptr+i) ).react(step, con_rates, error);
	}

	int new_count = 0; // new copy number `slash' new system array indicator

	for(int i=0; i<copyNum; ++i){
		switch ( react_vector[i] ) { // based on the reaction indicator create the system array
			case 0: // no reaction
				*(system_ptr+new_count) = old_system[i];
				new_count += 1;
				break;
			case 1: // replication
				*(system_ptr+new_count) = agent(current_id, (old_system[i]).get_id(), (old_system[i]).get_type(), (old_system[i]).get_rates());
				
				*(system_ptr+new_count+1) = agent(current_id+1, (old_system[i]).get_id(), (old_system[i]).get_type(), (old_system[i]).get_rates());
				new_count += 2;
				current_id += 2;
				break;
			case 3: // mutation
				*(system_ptr+new_count) = agent(current_id, (old_system[i]).get_id(), 0, wld_rates); // birth wild
				*(system_ptr+new_count+1) = agent(current_id+1, (old_system[i]).get_id(), 1, mnt_rates); // birth mutant
				new_count += 2;
				current_id += 2;
				break;
			default:
				break;
		}
	}
	for(int i=new_count; i<MAX_LEN; ++i)
		*(system_ptr+i) = agent();
}

void agented(int* output_array, agent* system_ptr, float Tmax, float step, float step_out, float init[2], float* wld_rates, float* mnt_rates, float* con_rates){

	int x[2];
	x[0] = init[0];
	x[1] = init[1];
	
	int copyNum = x[0]+x[1];
	int C0 = x[0] + x[1];
	int current_id = copyNum + 1; // (unique) id to be given to each new molecule in the system
	int Nout = (int) (Tmax/step_out + 1); // the number of time points saved in the output
	int Niter = (int) (Tmax/step + 1); // the number of iterations of the algorithm

	float target = step_out; // next time at which the output is saved
	float tt = 0; // current time
	int iter = 1; // interation number of the algorithm
	int out_iter = 1; // iteration of saved output
	
	*output_array = x[0];
	*(output_array+1) = x[1];
	
	for(int i=0; i<x[0]; ++i)
		*(system_ptr+i) = agent(i+1, 0, 0, wld_rates);
	for(int i=x[0]; i<copyNum; ++i)
		*(system_ptr+i) = agent(i+1, 0, 1, mnt_rates);

	while( tt<Tmax ){ // simulate forward in time
		agent old_state[copyNum];
		for(int i=0; i<copyNum; ++i)
			old_state[i] = *(system_ptr + i);
		
		int error = copyNum - C0;
		
		system_update(system_ptr, copyNum, error, step, current_id, wld_rates, mnt_rates, con_rates);
		
		bool mol_exist = true;
		x[0] = 0;
		x[1] = 0;
		int count = 0;
		
		while( mol_exist ){
			switch( (*(system_ptr+count)).get_type() ){
				case 0:
					x[0] += 1;
					count += 1;
					break;
				case 1:
					x[1] += 1;
					count += 1;
					break;
				default:
					mol_exist = false;
			}
		}
		copyNum = x[0] + x[1];
		assert( copyNum<=MAX_LEN);
		if( tt>target ){
			assert( out_iter<=Nout );
			*(output_array+2*out_iter) = x[0];
			*(output_array+2*out_iter+1) = x[1];
			out_iter += 1;
			target += step_out;
		}
		tt += step;
		iter += 1;
	}
}

int main(){
	filesystem::current_path( "/Users/jordanchilds/Documents/C++/ABM" );
	// global variables I like
	const float hour = 3600;
	const float day = 24*hour;
	const float year = 365*day;
	
	const int n_species = 2;
	
	float wild_rates[3] = {3.06e-8, 3.06e-8, 0};
	float mutant_rates[3] = {3.06e-8, 3.06e-8, 0.0};
	float control_rates[2] = {2.0e-3, 2.0e-3};
	
	float tmax = 100*year;
	float delta_t = 10*day;
	float stepOut = 1*year;
	int Nout = (int) (tmax/stepOut);
	float x_init[2] = {500,500};
	
	int output_array[Nout][n_species];
	int* output_ptr = &output_array[0][0];
	agent system_state[MAX_LEN];
	agent* system_ptr = &system_state[0];
	
	agented(output_ptr, system_ptr, tmax, delta_t, stepOut, x_init, wild_rates, mutant_rates, control_rates);
	
	for(int i=0; i<Nout; ++i){
		cout<< output_array[i][0] << " " << output_array[i][0] << endl;
	}
}

