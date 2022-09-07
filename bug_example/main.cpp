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

std::vector<Program> buildGraphAndPrograms( poplar::Graph &graph, long unsigned int x_len, const int numberOfCores, const int numberOfTiles, const int threadsPerTile) {
	
	long unsigned int datasetSize = numberOfCores*numberOfTiles*threadsPerTile ;
	int tileInt;

	// SHOULD WE PRE-COMPILE GILLESPIED? HOW YOU DO THAT?
	graph.addCodelets("example_codelet.cpp");
	
	Tensor vec = graph.addVariable(FLOAT, {x_len}, "a");
	
	Tensor output = graph.addVariable(FLOAT, {datasetSize, x_len}, "output");
	
	ComputeSet computeSet = graph.addComputeSet("computeSet");
	
	// Map tensors to tiles
	for(int i=0; i<datasetSize; ++i){
		
		int roundCount = i % int(numberOfCores * numberOfTiles * threadsPerTile);
		int tileInt = std::floor( float(roundCount) / float(threadsPerTile) );
		
		graph.setTileMapping(vec, tileInt);
		graph.setTileMapping(output[i], tileInt);
		
		VertexRef vtx = graph.addVertex(computeSet, "example_vertex");
		graph.setTileMapping(vtx, tileInt);
		
		graph.connect(vtx["x"], vec);
		graph.connect(vtx["out"], output[i]);
	}
	
	// to be able to read the output
	graph.createHostRead("output-read", output);
	
	// Create streams that allow reading and writing of the variables:
    //auto nTimes_stream = graph.addHostToDeviceFIFO("write_nTimes", INT, 1);
    auto stream = graph.addHostToDeviceFIFO("write_x", FLOAT, x_len);
	
	std::vector<Program> progs(Progs::NUM_PROGRAMS); // I HAVE NOT IDEA WHAT NUM_PROGRAMS IS/DOES - but without it the empire falls
	progs[WRITE_INPUTS] = Sequence( Copy(stream, vec.numElements()) );
	progs[CUSTOM_PROG] = Execute(computeSet);
						  
	return progs;
}

void executeGraphProgram(float* x_ptr, int x_len, poplar::Engine &engine) { // poplar::Device &device, std::vector<Program> progs, poplar::Graph &graph,

	engine.connectStream("write_x", x_ptr, x_ptr+x_len);
	engine.run(WRITE_INPUTS);
	engine.run(CUSTOM_PROG);
}


int main() {
	const int numberOfCores = 1; // access to POD16
	const int numberOfTiles = 1;// 1472;
	const int threadsPerTile = 1;
	
	long unsigned int totalThreads = numberOfCores*numberOfTiles*threadsPerTile ;
	
	auto manager = DeviceManager::createDeviceManager();
	auto devices = manager.getDevices(poplar::TargetType::IPU, numberOfCores);
	auto it = std::find_if(devices.begin(), devices.end(), [](Device &device) {
	return device.attach();
	});
	auto device = std::move(*it);
	
	Target target = device.getTarget();

	// Create the Graph object
	Graph graph(target);
	
	const int x_len = 5;
	float x[x_len] = {1,2,3,4,5};
	
	std::vector<Program> progs = buildGraphAndPrograms(graph, x_len, numberOfCores, numberOfTiles, threadsPerTile);
	
	Engine engine(graph, progs);
	engine.load(device);
	
	executeGraphProgram(x_ptr, x_len, engine);
						  
	std::vector<int> cpu_vector( totalThreads * x_len ) ;
	
	engine.readTensor("output-read", cpu_vector.data(), cpu_vector.data()+cpu_vector.size()) ;
	
	for(int i=0; i<x_len; ++i){
		cout << cpu_vector[ i ] << " ";
	}
	cout<<endl;

	return 0;
}
