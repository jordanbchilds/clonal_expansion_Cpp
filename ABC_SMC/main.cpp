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

double myMean(float* vec_ptr, long unsigned int len){
	double sum = std::accumulate(&vec_ptr[0], &vec_ptr[len], 0.0);
	double m =  sum / len;
}

double myMean(int* vec_ptr, long unsigned int len){
	double sum = std::accumulate(&vec_ptr[0], &vec_ptr[len], 0.0);
	double m =  (double) sum / (double) len;
}

double myStdDev(float* vec_ptr, long unsigned int len, float mean){
	double sq_diff = 0.0;
	for(int i=0; i<len; ++i)
		sq_diff += (*(vec_ptr + i)-mean) * (*(vec_ptr+i)-mean) ;

	return sq_diff / double(len - 1);
}

double myStdDev(int* vec_ptr, long unsigned int len, float mean){
	double sq_diff = 0.0;
	for(int i=0; i<len; ++i)
		sq_diff += (*(vec_ptr + i)-mean) * (*(vec_ptr+i)-mean) ;

	return sq_diff / double(len - 1);
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

std::vector<Program> buildGraphAndPrograms( poplar::Graph &graph, long unsigned int nParam, long unsigned int nTimes, const int numberOfCores, const int numberOfTiles, const int threadsPerTile) {
	
	long unsigned int totalThreads = numberOfCores*numberOfTiles*threadsPerTile ;
	int tileInt;

	// SHOULD WE PRE-COMPILE GILLESPIED? HOW YOU DO THAT?
	graph.addCodelets("gillespie_codelet.cpp");
	Tensor times = graph.addVariable(FLOAT, {nTimes}, "a");
	Tensor theta = graph.addVariable(FLOAT, {nParam}, "b");
	Tensor output = graph.addVariable(INT, {totalThreads, 2*nTimes}, "output");
	
	ComputeSet computeSet = graph.addComputeSet("computeSet");
	
	// Map tensors to tiles
	for(int i=0; i<totalThreads; ++i){
		int roundCount = i % int( totalThreads );
		int tileInt = std::floor( float(roundCount) / float(threadsPerTile) );
		
		graph.setTileMapping(times, tileInt);
		graph.setTileMapping(theta, tileInt);
		graph.setTileMapping(output[i], tileInt);
		
		VertexRef vtx = graph.addVertex(computeSet, "sim_network_vertex");
		graph.setTileMapping(vtx, tileInt);
		
		graph.connect(vtx["times"], times);
		graph.connect(vtx["theta"], theta);
		graph.connect(vtx["out"], output[i]);
	}
	
	// to be able to read the output
	graph.createHostRead("output-read", output);
	
	// Create streams that allow reading and writing of the variables:
    auto times_stream = graph.addHostToDeviceFIFO("write_dataTimes", FLOAT, nTimes);
	auto param_stream = graph.addHostToDeviceFIFO("write_theta", FLOAT, nParam);
	
	std::vector<Program> progs(Progs::NUM_PROGRAMS); // I HAVE NOT IDEA WHAT NUM_PROGRAMS IS/DOES - but without it the empire falls
	progs[WRITE_INPUTS] = Sequence( {Copy(param_stream, theta), Copy(times_stream, times)});
	progs[CUSTOM_PROG] = Execute(computeSet);
						  
	return progs;
}

void executeGraphProgram(float* theta_ptr, long unsigned int nParam, float* times_ptr, long unsigned int nTimes, poplar::Engine &engine) {

	engine.connectStream("write_dataTimes", times_ptr, times_ptr+nTimes);
	engine.connectStream("write_theta", theta_ptr, theta_ptr+nParam);
	
	engine.run(WRITE_INPUTS);
	engine.run(CUSTOM_PROG);
}


int main() {
	// READ IN THE DATA AND CALCULATE SUMMARY STATISTICS ARRAY
	// define input size - not ideal but we make do.

	 const long unsigned int nTimes = 3;
	 int nObs = 1000;
	
	 float ml_flat[nTimes*nObs];
	 int cn_flat[nTimes*nObs];
			
	 int count = 0;
	 fstream ml_file("./simulated_data/mutation_load.txt", ios_base::in);
	 if( ml_file.is_open() ){
		 float a;
		 while( ml_file >> a ){
			 ml_flat[count] = a;
			 count += 1;
		 }
	 }
	 ml_file.close();

	 count = 0;
	 fstream cn_file("./simulated_data/copy_number.txt", ios_base::in);
	 if( cn_file.is_open() ){
		 int a;
		 while( cn_file >> a ){
			 cn_flat[count] = a;
			 count += 1;
		 }
	 }
	 cn_file.close();
	 
	for(int i=0; i<nTimes*nObs; ++i){
		cout << cn_flat[i] << " " << ml_flat[i] << endl;
	}
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
		 data_summ[0][t][0] = myMean(&mut_load[t][0], nObs);
		 data_summ[0][t][1] = myStdDev(&mut_load[t][0], nObs, data_summ[0][t][0]);

		 data_summ[1][t][0] = myMean(&copy_num[t][0], nObs);
		 data_summ[1][t][1] = myStdDev(&copy_num[t][0], nObs, data_summ[1][t][0]);
	 }
	
	
	/*
	const int numberOfCores = 1; // access to POD16
	const int numberOfTiles = 1; // 1472;
	const int threadsPerTile = 1; // six threads per tile
	
	long unsigned int totalThreads = numberOfCores * numberOfTiles * threadsPerTile ;
	
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
	std::vector<Program> progs = buildGraphAndPrograms(graph, nParam, nTimes, numberOfCores, numberOfTiles, threadsPerTile);
	Engine engine(graph, progs);
	engine.load(device);

	
	 DEFINE PRIOR DISTRIBUTIONS
	 only used for first sample
	 std::default_random_engine generator;
	 std::normal_distribution<float> rate_dist(2.64e-3,5e-4);
	 std::normal_distribution<float> mut_dist(2e-4, 5e-5);
	 std::normal_distribution<float> con_dist(2e-3, 5e-4);
	 std::uniform_real_distribution<float> ML_dist(0.2,0.6);
	 std::normal_distribution<float> CN_dist(1e3, 100);
	
	//float param_space[Ntheta][nParam];
	
	GENERATE INITIAL PROPOSED PARAMETERS
	fill param_space array with initial params
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
	
	FOR ABC SMC (one day you'll get there, stay positive!)
	const int Nabc = 10;
	float thresholds[Nabc];
	float weights[Ntheta];
	myTheta param_state[Ntheta];
	for(int i=0; i<Nabc; ++i){
	weights[i] = 1.0;
	}
	
	//float threshold;
	
	
	float times[nTimes] = {25.0*365.0, 55.0*365.0, 65.0*365.0};
	float theta[nParam]  = {500.0, 500.0, 2.64e-3, 2.64e-3, 2.64e-3, 2.64e-3, 0.0, 2e-3, 2e-3};
	float* theta_ptr = &theta[0];
	float* times_ptr = &times[0];
	int Ntheta = 1;

	for(int i=0; i<Ntheta; ++i){
		for(int k=0; k<nParam; ++k){
			*(theta_ptr+k) = theta[k] ; //param_space[i][k];
		}
	
		executeGraphProgram(theta_ptr, nParam, times_ptr, nTimes, engine);
		
		std::vector<int> cpu_vector( totalThreads * 2*nTimes ) ;
		engine.readTensor("output-read", cpu_vector.data(), cpu_vector.data()+cpu_vector.size()) ;
		
		float mutation_load[nTimes][totalThreads];
		float copy_number[nTimes][totalThreads];
		float sim_summ[2][nTimes][2];
		for(int t=0; t<nTimes; ++t){
			for(int j=0; j<totalThreads; ++j){
				copy_number[t][j] = cpu_vector[j*2*nTimes + 2*t] + cpu_vector[j*2*nTimes + 2*t + 1];
				mutation_load[t][j] = cpu_vector[j*2*nTimes + 2*t + 1] / copy_number[t][j];
			}
			sim_summ[0][t][0] = myMean(&mutation_load[t][0], totalThreads);
			sim_summ[0][t][1] = myStdDev(&mutation_load[t][0], totalThreads, sim_summ[0][t][0]);
			sim_summ[1][t][0] = myMean(&copy_number[t][0], totalThreads);
			sim_summ[1][t][1] = myStdDev(&copy_number[t][0], totalThreads, sim_summ[1][t][0]);
		}
		double d = squared_dist(&sim_summ[0][0][0], &data_summ[0][0][0], nTimes, 2);
		cout<< d << endl;
	}
	*/
	return 0;
}
