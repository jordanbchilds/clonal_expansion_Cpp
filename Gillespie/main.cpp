//
//  main.cpp
//  Stochastic Simulation
//
//  Created by JORDAN CHILDS on 26/07/2022.
//
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>

using namespace std;

// global variables for mtDNA model
const int n_species = 2;
const int n_reactions = 5;
int Pre[n_reactions][n_species] = { {1,0}, {0,1}, {1,0}, {0,1}, {1,0} };
int Post[n_reactions][n_species] = { {2,0}, {0,2}, {0,0}, {0,0}, {1,1} };
float S[n_reactions][n_species] = { 0 };


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

double rand_unif(float lower=0, float upper=1){
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

float rep_controller(float con_rates[2], float rep_rate, float error) {
	if( error >= 0){
		return rep_rate*2.0/(1.0+exp(error*con_rates[0]));
	} else {
		return rep_rate*(1.0+exp(-1*error*con_rates[1]))/2.0;
	}
}

float* myHazards(float x[n_species], float rates[n_reactions], float error, float con_rates[2]){
	float temp_rates[n_reactions];
	temp_rates[0] = rep_controller(con_rates, rates[0], error);
	temp_rates[1] = rep_controller(con_rates, rates[1], error);
	for(int i=2; i<n_reactions; ++i)
		temp_rates[i] = rates[i];
	
    static float hazards[n_reactions];
    for(int i=0; i<n_reactions; ++i){
        float h_i = temp_rates[i];
        for(int j=0; j<n_species; ++j)
            h_i *= choose(x[j], Pre[i][j]);
        hazards[i] = h_i;
    }
    return hazards;
}

float* gillespied(float Tmax, float step_out, float x[n_species], float rates[n_reactions], float con_rates[2], float S[n_reactions][n_species]){
	// could use pointers and such instead of stating array size in function argument?
    // Tmax, step_out, x, rates all non-negative - use unsigned?
	int Nout = (int) (Tmax/step_out + 1);
    
    float out_arr[Nout][n_species];
    for(int j=0; j<n_species; ++j)
        out_arr[0][j] = x[j];
    
    // count, target, tt, are all non-negative - use unsigned?
    int count = 1;
    float target = step_out;
    float tt = 0;
    float* haz_ptr;
	float C0 = x[0]+x[1];
    
    while( tt<Tmax ){
		float copyNum = x[0]+x[1];
        haz_ptr = myHazards(x, rates, copyNum-C0, con_rates);
        float haz_total = 0; // hazard total non-negative - unsigned?
		
        for(int i=0; i<n_reactions; ++i)
            haz_total += *(haz_ptr + i);
        
        if( haz_total<1e-10 )
            tt = 1e10;
        else
            tt += rand_exp(haz_total);
        
        while( tt>=target ){
			//cout<< x[0]<<" "<<x[1] <<endl;
            for(int j=0; j<n_species; ++j)
                out_arr[count][j] = x[j];
            count += 1;
            target += step_out;
            if( count>Nout )
                goto save_output;
        }
        int r = rand_disc(n_reactions, haz_ptr);
		
        for(int i=0; i<n_species; ++i)
			x[i] += S[r][i];
    }
    save_output:
    /*
     ofstream myfile;
     myfile.open("Gillespie/Output/simmy_1.txt", ios::out);
     for(int row=0; row<Nout; row++){
         for(int col=0; col<n_species; col++)
             myfile << *(out_arr[row] + col)<<  " " ;
         myfile << "\n";
     }
     myfile.close();
     */
    return x;
}

int main() {
	// global variables I like
	float hour = 3600;
	float day = 24*hour;
	float year = 365*day;
    filesystem::current_path( "/Users/jordanchilds/Documents/C++/Gillespie" );
    
	for(int i=0; i<n_reactions; ++i){
		for(int j=0; j<n_species; ++j){
			S[i][j] = Post[i][j] - Pre[i][j];
		}
	}

    float x_init[n_species] = {500,500};
    float tmax = 100*year;
    float stepOut = 1*year;
    float react_rates[n_reactions] = { 3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0.0};
	float con_rates[2] = {0.0, 0.0};
	
	float* xx;
	srand48((long)time(NULL));
	xx = gillespied(tmax, tmax, x_init, react_rates, con_rates, S);
	cout << *xx << " " << *(xx+1) << "\n" ;
	
	for(int ii=0; ii<100; ++ii){
		float* xx;
		srand48((long) 101010101);
		xx = gillespied(tmax, tmax, x_init, react_rates, con_rates, S);
		cout << *xx << " " << *(xx+1) << "\n" ;
	}
	/*
	std::ofstream outfile ("test.txt");
	for(int ii=0; ii<1000; ++ii){
		float* xx;
        xx = gillespied(tmax, tmax, x_init, react_rates, con_rates, S);
		outfile << *xx << " " << *(xx+1) << "\n" ;
	}
	outfile.close();
	*/
}
