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

#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>
#include <poplar/DeviceManager.hpp>
#include <chrono>
#include <ctime> 

using namespace poplar;
using namespace poplar::program;


int main()
{
	const int numberOfTiles = 1472;
	const int threadsPerTile = 6;

	long unsigned int datasetSize = 8832; // anything less than 8832 takes the same runtime, 8833 takes double the runtime
	int tileInt;
	int nVal = 10000; // number of simulations

	float tmax = 3784320000.0; //120*365*24*60*60 in seconds
	float stepOut = 365.0*24.0*60.0*60.0; // in seconds
	long unsigned int Nout = (int) (tmax/stepOut + 1.0);
	
	// Create the DeviceManager which is used to discover devices
	auto manager = DeviceManager::createDeviceManager();
	// Attempt to attach to a single IPU:
	auto devices = manager.getDevices(poplar::TargetType::IPU, 1);
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
	 the input matrices of initial populations and rates have been split into vectors of length
	 datasetSize.
	 
	 Having them as matrices kept throwing errors that were making me upset.
	 */

	int w_initVals[datasetSize];
	int m_initVals[datasetSize];
	float reactOne_ratesVals[datasetSize];
	float reactTwo_ratesVals[datasetSize];
	float reactThree_ratesVals[datasetSize];
	float reactFour_ratesVals[datasetSize];
	float reactFive_ratesVals[datasetSize];
	float conOne_ratesVals[datasetSize];
	float conTwo_ratesVals[datasetSize];
  
	int w_initBoss = 500;
	int m_initBoss = 500;
    float reactOne_ratesBoss = 3.06e-8;
	float reactTwo_ratesBoss = 3.06e-8;
	float reactThree_ratesBoss = 3.06e-8;
	float reactFour_ratesBoss = 3.06e-8;
	float reactFive_ratesBoss = 0.0;
	float conOne_ratesBoss = 2.0e-3;
	float conTwo_ratesBoss = 2.0e-3;
  
    for (int i = 0; i < datasetSize; ++i){
		w_initVals[i] = w_initBoss;
		m_initVals[i] = m_initBoss;
		reactOne_ratesVals[i] =  reactOne_ratesBoss;
		reactTwo_ratesVals[i] =  reactTwo_ratesBoss;
		reactThree_ratesVals[i] =  reactThree_ratesBoss;
		reactFour_ratesVals[i] =  reactFour_ratesBoss;
		reactFive_ratesVals[i] =  reactFive_ratesBoss;
		conOne_ratesVals[i] = conOne_ratesBoss;
		conTwo_ratesVals[i] = conTwo_ratesBoss;
    }
	
	// Add steps to initialize the variables
	Tensor w_init = graph.addConstant<int>(INT, {datasetSize}, w_initVals);
	Tensor m_init = graph.addConstant<int>(INT, {datasetSize}, m_initVals);

	Tensor reactOne_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactOne_ratesVals);
	Tensor reactTwo_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactTwo_ratesVals);
	Tensor reactThree_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactThree_ratesVals);
	Tensor reactFour_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactFour_ratesVals);
	Tensor reactFive_rates = graph.addConstant<float>(FLOAT, {datasetSize}, reactFive_ratesVals);

	Tensor conOne_rates= graph.addConstant<float>(FLOAT, {datasetSize}, conOne_ratesVals);
	Tensor conTwo_rates= graph.addConstant<float>(FLOAT, {datasetSize}, conTwo_ratesVals);
	
	Tensor w_popDyn = graph.addVariable(INT, {datasetSize,Nout}, "w_popDyn");
	Tensor m_popDyn = graph.addVariable(INT, {datasetSize,Nout}, "m_popDyn");
	
	Tensor output = graph.addVariable(FLOAT, {datasetSize}, "output");

	ComputeSet computeSet = graph.addComputeSet("computeSet");

	// iterate through tiles on the IPU, map simulations to each thread (6) on each tile
	for (int i = 0; i < datasetSize; ++i)
	{
		int roundCount = i % int(numberOfTiles * threadsPerTile);
		int tileInt = std::floor(
			float(roundCount) / float(threadsPerTile)
			);
		// std::cout << "theadsPerIPU=" << int(numberOfTiles * threadsPerTile) << "i="<< i << " roundCount=" << roundCount << " tileInt=" << tileInt << std::endl;

		graph.setTileMapping(w_init[i], tileInt);
		graph.setTileMapping(m_init[i], tileInt);
		graph.setTileMapping(reactOne_rates[i], tileInt);
		graph.setTileMapping(reactTwo_rates[i], tileInt);
		graph.setTileMapping(reactThree_rates[i], tileInt);
		graph.setTileMapping(reactFour_rates[i], tileInt);
		graph.setTileMapping(reactFive_rates[i], tileInt);
		
		graph.setTileMapping(conOne_rates[i], tileInt);
		graph.setTileMapping(conTwo_rates[i], tileInt);
		
		graph.setTileMapping(w_popDyn[i], tileInt);
		graph.setTileMapping(m_popDyn[i], tileInt);
		
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
		graph.connect(vtx["w_popDyn"], w_popDyn[i]);
		graph.connect(vtx["m_popDyn"], m_popDyn[i]);
		graph.connect(vtx["out"], output[i]);

	}

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

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "Completed computation at " << std::ctime(&end_time)
			<< "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

	std::cout << "Rate=" << (float(datasetSize)/elapsed_seconds.count()) << std::endl;
	
	for(int i=0; i<datasetSize; ++i){
		cout<< output[i] <<endl;
	}
	/*
	ofstream myfile.open("ipu_copyNum.txt");
	for(int i=0; i<datasetSize; ++i)
		myfile<<output[i]<<endl;
	myfile.close();
	*/
	return 0;
}
