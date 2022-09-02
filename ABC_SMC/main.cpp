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
#include <boost/math/distributions/beta.hpp>

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
myTheta check_limits(myTheta theta, float* limits_array){
	bool inSupport = true;
	for(int i=0; i<9; ++i){
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
	theta_ss.muation = theta_star.mutation ;
	
	theta_ss.con_above = theta_star.con_above + runif(1e-4);
	theta_ss.con_abouse = theta_star.con_below +runif(1e-4);
	
	theta_ss.wInit = theta_star.wInit + runif_disc(25) ;
	theta_ss.mInit = theta_star.wInit + runif_disc(25) ;
	
	return theta_ss;
}

void create_graph(float* theta, long unsigned datasetSize){
	// iterate through tiles on the IPU, map simulations to each thread (6) on each tile
	for (int i = 0; i < datasetSize; ++i)
	{
		int roundCount = i % int(numberOfCores * numberOfTiles * threadsPerTile);
		int tileInt = std::floor( float(roundCount) / float(threadsPerTile) );
		graph.setTileMapping(w_init[i], tileInt);
		graph.setTileMapping(m_init[i], tileInt);
		graph.setTileMapping(reactOne_rates[i], tileInt);
		graph.setTileMapping(reactTwo_rates[i], tileInt);
		graph.setTileMapping(reactThree_rates[i], tileInt);
		graph.setTileMapping(reactFour_rates[i], tileInt);
		graph.setTileMapping(reactFive_rates[i], tileInt);
		
		graph.setTileMapping(conOne_rates[i], tileInt);
		graph.setTileMapping(conTwo_rates[i], tileInt);

		graph.setTileMapping(output[i], tileInt);

		VertexRef vtx = graph.addVertex(computeSet, "sim_network_vertex");
		graph.setTileMapping(vtx, tileInt);

		graph.connect(vtx["w_init"], w_init[i]);
		graph.connect(vtx["m_init"], m_init[i]);
		graph.connect(vtx["reactOne_rates"], reactOne_rates[i]);
		graph.connect(vtx["reactTwo_rates"], reactTwo_rates[i]);
		graph.connect(vtx["reactThree_rates"], reactThree_rates[i]);
		graph.connect(vtx["reactFour_rates"], reactFour_rates[i]);
		graph.connect(vtx["reactFive_rates"], reactFive_rates[i]);
		graph.connect(vtx["conOne_rates"], conOne_rates[i]);
		graph.connect(vtx["conTwo_rates"], conTwo_rates[i]);

		graph.connect(vtx["out"], output[i]);
	}

	// to be able to read the output
	graph.createHostRead("output-read", output);
	
	// Add a step to execute the compute set
	prog.add(Execute(computeSet));
	// Add a step to print out sim results
	// prog.add(PrintTensor("output", output));
	// Create the engine
	Engine engine(graph, prog);
	engine.load(device);

	auto start = std::chrono::system_clock::now();
	// Run the control program
	engine.run(0);
	auto end = std::chrono::system_clock::now();
}

int main()
{
	const int numberOfCores = 16; // access to POD16
	const int numberOfTiles = 1472;
	const int threadsPerTile = 6;

	 long unsigned int datasetSize = numberOfCores*numberOfTiles*threadsPerTile ; // 16 cores with 1472 tiles with 6 threads = 141,321 simulataneous simulations (thats a whole lotta simulations)
	 
	 int tileInt;

	 float tmax = 120.0*365.0*24.0*3600.0; // 120 years in seconds
	 float stepOut = 365.0*24.0*3600.0; // 1 year in seconds
	 long unsigned int Nout = (int) (tmax/stepOut + 1.0);
	 
	 // Create the DeviceManager which is used to discover devices
	 auto manager = DeviceManager::createDeviceManager();
	 // Attempt to attach to a single IPU:
	 auto devices = manager.getDevices(poplar::TargetType::IPU, numberOfCores);
	 
	 std::cout << "Trying to attach to IPU\n";
	 auto it = std::find_if(devices.begin(), devices.end(), [](Device &device) {
		 return device.attach();
	});
	 if (it == devices.end()) {
		 std::cerr << "Error attaching to device\n";
		 return -1;
	 }
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
	/*
	Add steps to initialize the variables
	Tensor w_init = graph.addVariable(INT, {datasetSize}, w_initVals);
	Tensor m_init = graph.addVariable(INT, {datasetSize}, m_initVals);
	Tensor reactOne_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactOne_ratesVals);
	Tensor reactTwo_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactTwo_ratesVals);
	Tensor reactThree_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactThree_ratesVals);
	Tensor reactFour_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactFour_ratesVals);
	Tensor reactFive_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactFive_ratesVals);
	
	Tensor conOne_rates= graph.addConstant<float>(FLOAT, {datasetSize}, conOne_ratesVals);
	Tensor conTwo_rates= graph.addConstant<float>(FLOAT, {datasetSize}, conTwo_ratesVals);
	 */
	Tensor popDyn = graph.addVariable(INT, {datasetSize, 2}, "popDyn") ;
	Tensor react_rates = graph.addVariable(FLOAT, {datasetSize, 5,2}, "react_rates") ;
	Tensor con_rates = graph.addVariable(FLOAT, {datasetSize, 2}, "con_rates");
	
	Tensor output = graph.addVariable(INT, {datasetSize, 2*Nout}, "output");
	/*
	ComputeSet computeSet = graph.addComputeSet("computeSet");
	const int Nabc = 10;
	float thresholds[Nabc];
	float weights[Ntheta];
	myTheta param_state[Ntheta];

	for(int i=0; i<Nabc; ++i){
		weights[i] = 1.0;
		// sample from some priors for inital state
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
