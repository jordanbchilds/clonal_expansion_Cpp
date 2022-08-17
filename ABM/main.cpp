#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// global variables I like
const float hour = 3600;
const float day = 24*hour;
const float year = 365*day;

// global variables for mtDNA model
const int MAX_LEN = 2000;
const int n_species = 2;

float rand_unif(float lower=0, float upper=1){
	srand( (unsigned)time(NULL) );
	rand();
	float unif_01;
	float unif;
	unif_01 = (float) rand() / (float) RAND_MAX;
	unif = unif_01*(upper-lower) + lower;
	return unif;
}

/*
#include <stdlib.h>double randZeroToOne(){
	return rand() / (RAND_MAX + 1.);
}
 */

class agent {
private:
	int unique_id;
	int parent_id;
	int species; // wild type: 0, mutant type: 1 (just like our matrix.)
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
			n_reactions = 2;
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

// constructor ??
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
	if(mol_type==0){
		n_reactions = 3;
	} else {
		n_reactions = 2;
	}
}
// deconstructor ??
agent::~agent(){}

/*
float rep_controller(float con_rates[2], float rep_rate, float error) {
	if( error>= 0){
		float out = rep_rate*2.0/(1.0+exp(error*con_rates[0]));
		return out;
	} else {
		float out = rep_rate*(1.0+exp(-1*error*con_rates[1]))/2.0;
		return out;
	}
}
 */
float rep_controller(float con_rates[2], float rep_rate, float error){
	return 3.06e-8;
}

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
		haz_total += *(rates+i); // rates should be a pointer to a vector of length n_reactions
	if( haz_total==0 )
		return 0;
	
	float event_prob = 1 - exp(-1*haz_total*deltaT);
	//double u1 = rand_unif();
	float u1 = drand48();
	/*
	 the probability of an event occurring in (t,t+dt] is given as 1-exp(-h0*dt), where h0 is the total hazard
	 */
	if( u1<event_prob ){
		// could change this to use the rand_disc function you wrote? - in Gillespie file
		float cumProb[n_reactions];
		cumProb[0] = *rates /haz_total;
		
		for(int i=1; i<n_reactions; ++i){
			float p_i = *rates /haz_total;
			for(int j=1; j<=i; ++j)
				p_i += *(rates+j)/haz_total;
			cumProb[i] = p_i;
		}
		
		//double u2 = rand_unif();
		float u2 = drand48();
		if( u2<=cumProb[0] )
			return 1; // replication: 1
		for(int i=1; i<n_reactions; ++i){
			if(cumProb[i-1]<u2 && u2<=cumProb[i])
				return i+1; // degradation: 2, mutation: 3
		}
	}
	return 0; // no reaction: 0
}

agent* system_update(agent* current_ptr, float step, int current_id, float* wld_rates, float* mnt_rates, float* con_rates){
	
	int react_vector[MAX_LEN]; // array to store reaction that occur for each mol 0,1,...,(MAX_LEN-1)

	// pararallellelise
	for( int i=0; i<MAX_LEN; ++i)
		react_vector[i] = (*(current_ptr+i)).react(step, con_rates);
	
	static agent new_system[MAX_LEN]; // saves to stack so can be accessed outside function
	bool mol_exist = true; // boolean to indicate if the molecule in the massive array exists or if it is just a place holder
	
	int new_count = 0; // new copy number `slash' new system array indicator
	int old_count = 0; // old copy number `slash' old system array indicator
	
	while( mol_exist ){ // loop to create to system array
		switch ( react_vector[old_count] ) { // based on the reaction indicator create the system array
			case 0: // no reaction
				new_system[new_count] = *(current_ptr+old_count);
				new_count += 1;
				old_count += 1;
				break;
			case 1: // replication
				new_system[new_count] = agent(current_id, (*(current_ptr+old_count)).get_id(), (*(current_ptr+old_count)).get_type(), (*(current_ptr+old_count)).get_rates());
				
				new_system[new_count+1] = agent(current_id+1, (*(current_ptr+old_count)).get_id(), (*(current_ptr+old_count)).get_type(), (*(current_ptr+old_count)).get_rates());
				old_count += 1;
				new_count += 2;
				current_id += 2;
				break;
			case 2: // degradation
				old_count += 1;
				break;
			case 3: // mutation
				new_system[new_count] = agent(current_id, (*(current_ptr+old_count)).get_id(), 0, wld_rates); // birth wild
				new_system[new_count+1] = agent(current_id+1, (*(current_ptr+old_count)).get_id(), 1, mnt_rates); // birth mutant
				old_count += 1;
				new_count += 2;
				current_id += 2;
				break;
			default:
				// if the molecule is purely a place holder in our massive system array the reat() output should be one
				// when the first react() = -1 occurs we can exist our loop
				mol_exist = false;
		}
	}
	return new_system; // return new system state
}

float* agented(float Tmax, float step, float step_out, float init[n_species], float* wld_rates, float* mnt_rates, float* con_rates){
	// define population vector, x to be updated with each iteration
	int x[n_species];
	x[0] = init[0];
	x[1] = init[1];
	
	int copyNum = x[0] + x[1];
	int current_id = copyNum + 1; // (unique) id to be given to each new molecule in the system
	int Nout = (int) (Tmax/step_out + 1); // the number of time points saved in the output
	int Niter = (int) (Tmax/step + 1); // the number of iterations of the algorithm

	float out_arr[Nout][2]; // array to store the output
	float target = step_out; // next time at which the output is saved
	float tt = 0; // current time
	int iter = 1; // interation number of the algorithm
	int out_iter = 1;
	
	out_arr[0][0] = x[0]; // input the starting conditions to the output array
	out_arr[0][1] = x[1]; // ditto
	
	agent system_state[MAX_LEN]; // create a system array of length MAX_LEN
	
	// create the molecules and save them in the system array
	for(int i=0; i<x[0]; ++i)
		*(system_state+i) = agent(i+1, 0, 0, wld_rates);
	for(int i=x[0]; i<copyNum; ++i)
		*(system_state+i) = agent(i+1, 0, 1, mnt_rates);

	while( tt<Tmax ){ // simulate forward in time
		agent old_state[MAX_LEN]; // define
		for(int i=0; i<copyNum; ++i)
			*(old_state + i) = *(system_state + i);
			
		agent* system_state;
		system_state = system_update(old_state, step, current_id, wld_rates, mnt_rates, con_rates);
		
		bool mol_exist = true;
		int copyNum = 0;
		x[0] = 0;
		x[1] = 0;
		while( mol_exist ){
			switch ((*(system_state+copyNum)).get_type()) {
				case 1:
					x[0] += 1;
					copyNum +=1 ;
					break;
				case 0:
					x[1] += 1;
					copyNum += 1;
					break;
				default:
					mol_exist = false;
			}
		}

		while( tt>target ){
			assert( out_iter<=Nout );
			// cout<< x[0] << " " <<x[1] <<endl;
			out_arr[out_iter][0] = x[0];
			out_arr[out_iter][1] = x[1];
			out_iter += 1;
			target += step_out;
		}
		tt += step;
		iter += 1;
		assert( iter<=Niter );
	}

}

int main(){
	srand48((long)time(NULL));
	
	float wild_rates[3] = {3.06e-8, 3.06e-8, 3.06e-11};
	float mutant_rates[2] = {3.06e-5, 3.06e-5};
	float control_rates[2] = {5e-4, 5e-4};
	
	float tmax = 100*year;
	float delta_t = 10*day;
	float stepOut = 50*day;
	float x_init[2] = {100,100};
	
	for(int i=0; i<100; ++i){
		agented(tmax, delta_t, stepOut, x_init, wild_rates, mutant_rates, control_rates);
	}

	
}

