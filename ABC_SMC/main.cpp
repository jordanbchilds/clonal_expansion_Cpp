// Copyright (c) 2022 Graphcore Ltd. All rights reserved.
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>
#include <poplar/DeviceManager.hpp>
#include <chrono>
#include <ctime>

using namespace std;
using namespace poplar;
using namespace poplar::program;

int runif_disc(unsigned h){
	return  rand() % h;
}

float runif_01(){
	return (float) rand() / (float) RAND_MAX ;
}

float runif(float h){
	float unif_01 = (float) rand() / (float) RAND_MAX;
	return 2.0*unif_01*h - h ;
}

float runif(float lower, float upper){
	float unif_01 = (float) rand() / (float) RAND_MAX;
	return unif_01*(upper - lower) + lower ;
}

int rdisc(unsigned n, float* weights){
	float norm = 0.0; // normalising constant
	float cumWeights[n];
	for(int i=0; i<n; ++i)
		norm += *(weights+i);
	
	for(int i=0; i<n; ++i){
		float cc = 0.0;
		for(int j=0; j<=i; ++j)
			cc += *(weights+j)/norm;
		cumWeights[i] = cc;
	}

	float u = runif_01() ;
	if( 0<=u && u<cumWeights[0] ){
		return 0;
	} else {
		for(int i=1; i<n; ++i){
			if( cumWeights[i-1]<u && u<=cumWeights[i] )
				return i;
		}
	}
	return -1; // requires a return outside the loop
}

/*
myTheta check_limits(myTheta theta, float* limits_array, int nParam){
	bool inSupport = true;
	for(int i=0; i<nParam; ++i){
		inSupport &= (
	}
	re
}
*/

void perturb( float* theta_star ){
	*(theta_star) = *theta_star + runif_disc(25) ;
	*(theta_star+1) = *(theta_star+1) + runif_disc(25) ;
	
	*(theta_star+2) = *(theta_star+2) + runif(1e-5);
	*(theta_star+3) = *(theta_star+3) + runif(1e-5);
	*(theta_star+4) = *(theta_star+4) + runif(1e-5);
	*(theta_star+5) = *(theta_star+5) + runif(1e-5);
	*(theta_star+6) = *(theta_star+6) + runif(1e-6) ;
	
	*(theta_star+7) = *(theta_star+7) + runif(1e-5);
	*(theta_star+8) = *(theta_star+8) + runif(1e-5);
}

double myMean(float* vec_ptr, int len){
	double sum = std::accumulate(&vec_ptr[0], &vec_ptr[len], 0.0);
	double m =  sum / len;
}

double myMean(int* vec_ptr, int len){
	double sum = std::accumulate(&vec_ptr[0], &vec_ptr[len], 0.0);
	double m =  sum / len;
}

double myStdDev(float* vec_ptr, int len, double mean){
		double accum = 0.0;
				std::for_each (vec_ptr, vec_ptr+len, [&](const double d) {
		   accum += (d - mean) * (d - mean);
	   });
	return sqrt(accum / len-1);
}

double squared_dist(float* x_ptr, float* y_ptr, int nrow, int ncol){
	
	double distance = 0.0;
	for(int i=0; i<nrow; ++i){
		for(int j=0; j<ncol; ++j){
			double xy_diff = *(x_ptr+i*ncol+j) - *(y_ptr+i*ncol+j);
			distance += xy_diff*xy_diff;
		}
	}
	return distance ;
}

enum Progs {
	WRITE_INPUTS,
	CUSTOM_PROG,
	NUM_PROGRAMS
};

std::vector<Program> buildGraphAndPrograms( poplar::Graph &graph, long unsigned int nParam, long unsigned int nTimes ) {
	
	const int numberOfCores = 16; // access to POD16
	const int numberOfTiles = 1472;
	const int threadsPerTile = 6;

	long unsigned int datasetSize = numberOfCores*numberOfTiles*threadsPerTile ;
	int tileInt;

	float tmax = 120.0*365.0; // 120 years in seconds
	float stepOut = 365.0; // 1 year in seconds

	// SHOULD WE PRE-COMPILE GILLESPIED? HOW YOU DO THAT?
	graph.addCodelets("gillespie_codelet.cpp");
	Tensor Nout = graph.addVariable(INT, {1}, "n_times");
	Tensor times = graph.addVariable(FLOAT, {nTimes}, "data_times");
	Tensor theta = graph.addVariable(FLOAT, {nParam}, "model_params");
	Tensor output = graph.addVariable(INT, {datasetSize, Nout*2}, "output");
	ComputeSet computeSet = graph.addComputeSet("computeSet");
	
	// Map tensors to tiles
	for(int i=0; i<datasetSize; ++i){
		int roundCount = i % int(numberOfCores * numberOfTiles * threadsPerTile);
		int tileInt = std::floor( float(roundCount) / float(threadsPerTile) );
		graph.setTileMapping(Nout, tileInt);
		graph.setTileMapping(times, tileInt);
		graph.setTileMapping(theta, tileInt);
		graph.setTileMapping(output[i], tileInt);
		
		VertexRef vtx = graph.addVertex(computeSet, "sim_network_vertex");
		
		graph.setTileMapping(vtx, tileInt);
		graph.connect(vtx["Nout"], Nout);
		graph.connect(vtx["times"], times);
		graph.connect(vtx["theta"], theta);
		graph.connect(vtx["out"], output[i]);
	}
	
	// to be able to read the output
	graph.createHostRead("output-read", output);
	
	// Create streams that allow reading and writing of the variables:
    auto nTimes_stream = graph.addHostToDeviceFIFO("write_nTimes", INT, 1);
    auto times_stream = graph.addHostToDeviceFIFO("write_dataTimes", FLOAT, nTimes);
	auto param_stream = graph.addHostToDeviceFIFO("write_theta", FLOAT, nParam);
	
	std::vector<Program> progs(Progs::NUM_PROGRAMS); // I HAVE NOT IDEA WHAT NUM_PROGRAMS IS/DOES - but without it the empire falls
	progs[WRITE_INPUTS] = Copy(param_stream, theta);
	progs[CUSTOM_PROG] = Execute(computeSet);
						  
	return progs;
}

void executeGraphProgram(float* theta_ptr, int nParam, float* outTimes_ptr, int nTimes, poplar::Engine &engine) { // poplar::Device &device, std::vector<Program> progs, poplar::Graph &graph,
	
	engine.connectStream("write_nTimes", &nTimes, &nTimes);
	engine.connectStream("write_dataTimes", outTimes_ptr, outTime_ptr+nTimes);
	engine.connectStream("write_theta", theta_ptr, theta_ptr+nParam);
	
	engine.run(WRITE_INPUTS);
	engine.run(CUSTOM_PROG);
}


int main() {
	// READ IN THE DATA AND CALCULATE SUMMARY STATISTICS ARRAY
	// define input size - not ideal but we make do.
	 long unsigned int nTimes = 3;
	 long unsigned int nObs = 1000;
	 
	 float ml_flat[nTimes*nObs];
	 int cn_flat[nTimes*nObs];
			
	 int count = 0;
	 fstream ml_file("./mutation_load.txt", ios_base::in);
	 if( ml_file.is_open() ){
		 float a;
		 while( ml_file >> a ){
			 ml_flat[count] = a;
			 count += 1;
		 }
	 }
	 ml_file.close();

	 count = 0;
	 fstream cn_file("./copy_number.txt", ios_base::in);
	 if( cn_file.is_open() ){
		 int a;
		 while( cn_file >> a ){
			 cn_flat[count] = a;
			 count += 1;
		 }
	 }
	 cn_file.close();
	 
	 int copy_num[nTimes][nObs];
	 float mut_load[nTimes][nObs];
	 for(int i=0; i<nTimes; ++i){
		for(int j=0; j<nObs; ++j){
			copy_num[i][j] = cn_flat[i*nObs +j];
			mut_load[i][j] = ml_flat[i*nObs +j];
		}
	 }

	 float data_summ[2][nTimes][2];
	 for(int t=0; t<nTimes; ++t){
		 data_summ[0][t][0] = myMean(mut_load[t], nObs);
		 data_summ[0][t][1] = myStdDev(mut_load[t], nObs, data_summ[0][t][0]);

		 data_summ[1][t][0] = myMean(copy_num[t], nObs);
		 data_summ[1][t][1] = myStdDev(copy_num[t], nObs, data_summ[1][t][0]);
	 }
	 
	const int numberOfCores = 16; // access to POD16
	const int numberOfTiles = 1472;
	const int threadsPerTile = 6;
	long unsigned int totalThreads = numberOfCores*numberOfTiles*threadsPerTile ;

	// const unsigned Nsim = 1e5;
	// int thetaIter = int(totalThreads/Nsim);

	auto manager = DeviceManager::createDeviceManager();
	auto devices = manager.getDevices(poplar::TargetType::IPU, numberOfCores);
	auto it = std::find_if(devices.begin(), devices.end(), [](Device &device) {
	return device.attach();
	});
	auto device = std::move(*it);

	std::cout << "Attached to IPU " << device.getId() << std::endl;
	Target target = device.getTarget();

	long unsigned int nParam = 9;
	
	// Create the Graph object
	Graph graph(target);
	std::vector<Program> progs = buildGraphAndPrograms(graph, nParam, nTimes);
	Engine engine(graph, progs);
	engine.load(device);
						  
	/*
	DEFINE PRIOR DISTRIBUTIONS
	only used for first sample
	*/
	std::default_random_engine generator;
	std::normal_distribution<float> rate_dist(2.64e-3,5e-4);
	std::normal_distribution<float> mut_dist(2e-4, 5e-5);
	std::normal_distribution<float> con_dist(2e-3, 5e-4);
	std::uniform_real_distribution<float> ML_dist(0.2,0.6);
	std::normal_distribution<float> CN_dist(1e3, 100);

	float times[nTimes] = {25.0,55.0,65.0};
	float Tmax = 120.0*365.0 ;
	float theta[nParam] ;
	float* theta_ptr = &theta[0];
	const unsigned Ntheta = 100;
						  
	float param_space[Ntheta][nParam];
	/*
	GENERATE INITIAL PROPOSED PARAMETERS
	fill param_space array with initial params
	*/
	for(int i=0; i<Ntheta; ++i){
		int c0 = CN_dist(generator);
		float h0 = ML_dist(generator);
		param_space[i][0] = round( c0*(1.0-h0) );
		param_space[i][1] = round( c0*h0 );

		param_space[i][2] = rate_dist(generator);
		param_space[i][3] = rate_dist(generator);
		param_space[i][4] = rate_dist(generator);
		param_space[i][5] = rate_dist(generator);
		param_space[i][6] = mut_dist(generator);
		param_space[i][7] = con_dist(generator);
		param_space[i][8] = con_dist(generator);
	}
	/* FOR ABC SMC (one day you'll get there, stay positive!)
	const int Nabc = 10;
	float thresholds[Nabc];
	float weights[Ntheta];
	myTheta param_state[Ntheta];
	for(int i=0; i<Nabc; ++i){
	weights[i] = 1.0;
	}
	*/
	// float threshold;

	for(int i=0; i<Ntheta; ++i){
		for(int k=0; k<nParam; ++k){ *(theta_ptr+k) = param_space[i][k]; }
				
		executeGraphProgram(theta_ptr, nParam, &times[0], nTimes, engine);
		std::vector<int> cpu_vector( totalThreads * nTimes * 2 );
		engine.readTensor("output-read", cpu_vector.data(), cpu_vector.data()+cpu_vector.size());

		float sim_summ[2][nTimes][2];
		for(int t=0; t<nTimes; ++t){
			float mutation_load[nTimes][totalThreads];
			float copy_number[nTimes][totalThreads];
			
			for(int i=0; i<totalThreads; ++i){
				copy_number[t][i] = cpu_vector[i*2*nTimes + 2*t] + cpu_vector[i*2*nTimes + 2*t + 1];
				mutation_load[t][i] = cpu_vector[i*2*nTimes + 2*t + 1] / copy_number[t][i];
			}
			sim_summ[0][t][0] = myMean(mutation_load[t], totalThreads);
			sim_summ[0][t][1] = myStdDev(mutation_load[t], totalThreads, sim_summ[0][t][0]);
			sim_summ[1][t][0] = myMean(copy_number[t], totalThreads);
			sim_summ[1][t][1] = myStdDev(copy_number[t], totalThreads, sim_summ[1][t][0]);
		}
		double d = squared_dist(sim_summ, data_summ, nTimes, 2);
		cout<< d <<endl;
	}
	
	return 0;
}
