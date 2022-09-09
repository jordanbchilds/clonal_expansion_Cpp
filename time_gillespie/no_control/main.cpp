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
	const int numberOfCores = 16; // access to POD16
	const int numberOfTiles = 1472; // 1472;
	const int threadsPerTile = 6; // six threads per tile
	
	long unsigned int totalThreads = numberOfCores * numberOfTiles * threadsPerTile ;
	
	auto manager = DeviceManager::createDeviceManager();
	auto devices = manager.getDevices(poplar::TargetType::IPU, numberOfCores);
	auto it = std::find_if(devices.begin(), devices.end(), [](Device &device) {
	return device.attach();
	});
	auto device = std::move(*it);

	std::cout << "Attached to IPU " << device.getId() << std::endl;
	Target target = device.getTarget();

	const long unsigned int nParam = 7;
	const long unsigned int nTimes = 100;
	
	// Create the Graph object
	Graph graph(target);
	std::vector<Program> progs = buildGraphAndPrograms(graph, nParam, nTimes, numberOfCores, numberOfTiles, threadsPerTile);
	Engine engine(graph, progs);
	engine.load(device);
	
	float times[nTimes+1];
	for(int i=0; i<=nTimes; ++i){
		times[i] = i*365.0*24.0*3600.0 ;
	}
	float theta[nParam] = {500.0, 500.0, 3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8}// = {500.0, 500.0, 2.64e-3, 2.64e-3, 2.64e-3, 2.64e-3, 0.0};
	float* theta_ptr = &theta[0];
	float* times_ptr = &times[0];
	
	const int Nsim = 10;
	double simTimes[Nsim] = {0};
	for(int t=0; t<Nsim; ++t){
		
		auto start = chrono::high_resolution_clock::now();
		executeGraphProgram(theta_ptr, nParam, times_ptr, nTimes, engine);
		auto end = chrono::high_resolution_clock::now();
		
		chrono::duration<double, std::milli> ms_double = end - start;
		
		simTimes[t] = ms_double.count() ;
		cout<< ms_double.count() << endl;
	}
	
	/*
	std::vector<int> cpu_vector( totalThreads * 2*nTimes ) ;
	engine.readTensor("output-read", cpu_vector.data(), cpu_vector.data()+cpu_vector.size()) ;
	*/
	
	return 0;
}
