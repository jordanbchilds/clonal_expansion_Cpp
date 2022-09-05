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

using namespace poplar;
using namespace poplar::program;

struct myTheta {
	float rep_wld;
	float deg_wld;
	float rep_mnt;
	float deg_mnt;
	float mutation;
	
	float con_above;
	float con_below;
	
	unsigned wInit;
	unsigned mInit;
} ;

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

myTheta perturb(myTheta theta_star){
	myTheta theta_ss;
	theta_ss.rep_wld = theta_star.rep_wld + runif(1e-4);
	theta_ss.rep_mnt = theta_star.rep_mnt + runif(1e-4);
	theta_ss.deg_wld = theta_star.deg_wld + runif(1e-4);
	theta_ss.deg_mnt = theta_star.deg_mnt + runif(1e-4);
	theta_ss.mutation = theta_star.mutation ;
	
	theta_ss.con_above = theta_star.con_above + runif(1e-4);
	theta_ss.con_below = theta_star.con_below +runif(1e-4);
	
	theta_ss.wInit = theta_star.wInit + runif_disc(25) ;
	theta_ss.mInit = theta_star.wInit + runif_disc(25) ;
	
	return theta_ss;
}

enum Progs {
	WRITE_INPUTS,
	CUSTOM_PROG,
	READ_RESULTS
};

poplar::program::Program buildGraphAndPrograms( poplar::Graph &graph ) {
	// Use the namespace here to make graph construction code less verbose:
	using namespace poplar;
	
	const int numberOfCores = 16; // access to POD16
	const int numberOfTiles = 1472;
	const int threadsPerTile = 6;

	long unsigned int datasetSize = numberOfCores*numberOfTiles*threadsPerTile ;
	int tileInt;

	float tmax = 120.0*365.0; // 120 years in seconds
	float stepOut = 365.0; // 1 year in seconds
	long unsigned int Nout = (int) (tmax/stepOut + 1.0);
	
	std::size_t nParam = 9;
	// Create Tensor of params to be passed to tile
	Tensor theta = graph.addVariable(FLOAT, {nParam}, "a");
	
	// Tensor to stoer output from each tile
	Tensor output = graph.addVariable(INT, {datasetSize, 2*Nout}, "output");
	// In order to do any computation we need a compute set and a compute
	// vertex that is placed in that compute set:
	ComputeSet computeSet = graph.addComputeSet("computeSet");
	
	// Map tensors to tiles
	for(int i=0; i<datasetSize; ++i){
		int roundCount = i % int(numberOfCores * numberOfTiles * threadsPerTile);
		int tileInt = std::floor( float(roundCount) / float(threadsPerTile) );
		
		graph.setTileMapping(theta, tileInt);
		
		VertexRef vtx = graph.addVertex(computeSet, "sim_network_vertex");
		graph.setTileMapping(vtx, tileInt);
		
		graph.connect(vtx["theta"], theta);
		graph.connect(vtx["out"], output[i]);
	}
		

	// SHOULD WE PRE-COMPILE GILLESPIED? HOW YOU DO THAT?
	graph.addCodelets("gillespie_codelets.cpp");

	// Create streams that allow reading and writing of the variables:
	auto param_stream = graph.addHostToDeviceFIFO("write_theta", FLOAT, nParam);
	auto output_inStream = graph.addHostToDeviceFIFO("write_output", FLOAT, output.numElements());
	auto output_outStream = graph.addDeviceToHostFIFO("read_output", FLOAT, output.numElements());
	// auto stream4 = g.addDeviceToHostFIFO("read_z",  FLOAT, v3.numElements());
	// I DON'T THINK I NEED AN OUTPUT STREAM - OUTPUT ALREADY OUTPUT'ING

	program::Program progs;

	// Add program which initialises the inputs. Poplar is able to merge these
	// copies for efficiency:
	progs[WRITE_INPUTS] = program::Sequence({program::Copy(input_stream, theta),program::Copy(output_inStream, output)});

	// Program that executes custom vertex in compute set 1:
	progs[CUSTOM_PROG] = program::Execute(computeSet);

	// Add a program to read back the result:
	progs[READ_RESULTS] = program::Copy(output, output_outStream);

	return progs;
}


int main()
{
	const int numberOfCores = 16; // access to POD16
	const int numberOfTiles = 1472;
	const int threadsPerTile = 6;

	long unsigned int datasetSize = numberOfCores*numberOfTiles*threadsPerTile ;
	auto manager = DeviceManager::createDeviceManager(); // Create the DeviceManager which is used to discover devices
	auto devices = manager.getDevices(poplar::TargetType::IPU, numberOfCores); // Attempt to attach to a single IPU
	auto it = std::find_if(devices.begin(), devices.end(), [](Device &device) {
		return device.attach();
	});
	auto device = std::move(*it);
	std::cout << "Attached to IPU " << device.getId() << std::endl;
	Target target = device.getTarget();

	// Create the Graph object
	Graph graph(target);

	// Add codelets to the graph
	graph.addCodelets("gillespie_codelet.cpp");

	// Create a control program that is a sequence of steps
	Sequence prog;
	
	/*
	DEFINE PRIOR DISTRIBUTIONS
	only used for first sample
	*/
	std::default_random_engine generator;
	std::normal_distribution<float> rate_dist(3.06e-8,5e-9);
	std::normal_distribution<float> con_dist(2e-3, 5e-4);
	std::uniform_real_distribution<float> ML_dist(0.45,0.55);
	std::normal_distribution<float> CN_dist(1e4, 250);

	const unsigned Ntheta = 1000;
	myTheta param_space[Ntheta];
	/*
	GENERATE INITIAL PROPOSED PARAMETERS
	fill param_space array with initial params
	*/
	for(int i=0; i<Ntheta; ++i){
		param_space[i].rep_wld = rate_dist(generator);
		param_space[i].rep_mnt = rate_dist(generator);
		param_space[i].deg_wld = rate_dist(generator);
		param_space[i].deg_mnt = rate_dist(generator);
		param_space[i].mutation = 0.0;

		param_space[i].con_above = con_dist(generator);
		param_space[i].con_below = con_dist(generator);

		int c0 = CN_dist(generator);
		float h0 = ML_dist(generator);
		param_space[i].wInit = round( c0*(1.0-h0) );
		param_space[i].mInit = round( c0*h0 );
	}
	
	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "Completed computation at " << std::ctime(&end_time)
	<< "Time to create "<< Ntheta<<" graphs: " << elapsed_seconds.count() << "s" << std::endl;
	
	/*
	ComputeSet computeSet = graph.addComputeSet("computeSet");
	const int Nabc = 10;
	float thresholds[Nabc];
	float weights[Ntheta];
	myTheta param_state[Ntheta];

	for(int i=0; i<Nabc; ++i){
	    weights[i] = 1.0;
	}

	for(int t=0; t<Nabc; ++t){
	    int i = 0;
		while( i<Ntheta ){
			unsigned index = rdisc(Ntheta, weights);
			myTheta theta_star = perturb(*(param_state+index));

		// create_graph();
		}
	}

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "Completed computation at " << std::ctime(&end_time)
	<< "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
	std::vector<int> cpu_vector( datasetSize * Nout * 2 );

	engine.readTensor("output-read", cpu_vector.data(), cpu_vector.data()+cpu_vector.size());

	std::ofstream wild_file ("ipu_wldCount.txt");
	for(int i=0; i<datasetSize; ++i){
		for(int j=0; j<Nout; ++j){
			wild_file<< cpu_vector[ i*2*Nout + 2*j ] << "\t" ;
		}
		wild_file<< "\n";
	}
	wild_file.close();

	std::ofstream mtnt_file ("ipu_mntCount.txt");
	for(int i=0; i<datasetSize; ++i){
		for(int j=0; j<Nout; ++j){
			mtnt_file<< cpu_vector[ i*2*Nout + 2*j + 1 ] << "\t" ;
		}
		mtnt_file<< "\n";
	}
	mtnt_file.close();
	*/
	return 0;
}
