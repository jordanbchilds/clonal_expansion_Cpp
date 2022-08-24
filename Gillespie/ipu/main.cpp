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

  int x_initVals[datasetSize][2];
  float react_ratesVals[datasetSize][5];
  float con_ratesVals[datasetSize][2];
  
  int x_initBoss[2] = {500,500};
  float react_ratesBoss[5] = { 3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0.0};
  float con_ratesBoss[2] = {2.0e-3, 2.0e-3};
  
  for (int i = 0; i < datasetSize; ++i){
	  for(int j=0; j<2; j++){
		x_initVals[i][j] = x_initBoss[j];
	  }
	  for(int j=0; j<5; j++){
		react_ratesVals[i][j] =  react_ratesBoss[j];
	  }
	  for(int j=0; j<2; j++){
		con_ratesVals[i][j] = con_ratesBoss[j];
	  }
  }

  // Add steps to initialize the variables
  Tensor x_init = graph.addConstant<int>(INT, {datasetSize}, x_initVals);
  Tensor react_rates = graph.addConstant<float>(FLOAT, {datasetSize,5}, react_ratesVals);
  Tensor con_rates= graph.addConstant<float>(FLOAT, {datasetSize,2}, con_ratesVals);

  Tensor out = graph.addVariable(INT, {datasetSize,2,int(tmax/step_out + 1.0)}, "output");

  ComputeSet computeSet = graph.addComputeSet("computeSet");

  // iterate through tiles on the IPU, map 6 sets of variables and 6 option pricers to each tile
  for (int i = 0; i < datasetSize; ++i)
  {
    int roundCount = i % int(numberOfTiles * threadsPerTile);
    int tileInt = std::floor(
        float(roundCount) / float(threadsPerTile)
        );
    // std::cout << "theadsPerIPU=" << int(numberOfTiles * threadsPerTile) << "i="<< i << " roundCount=" << roundCount << " tileInt=" << tileInt << std::endl;

    graph.setTileMapping(x_init[i], tileInt);
    graph.setTileMapping(react_rates[i], tileInt);
    graph.setTileMapping(con_rates[i], tileInt);
    graph.setTileMapping(out[i], tileInt);

    VertexRef vtx = graph.addVertex(computeSet, "sim_network");
    graph.setTileMapping(vtx, tileInt);

    graph.connect(vtx["x_init"], x_init[i]);
    graph.connect(vtx["react_rates"], react_rates[i]);
    graph.connect(vtx["con_rates"], con_rates[i]);
    graph.connect(vtx["out"], out[i]);

  }

  // Add a step to execute the compute set
  prog.add(Execute(computeSet));
  // Add a step to print out sim results
  prog.add(PrintTensor("out", out));
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

  return 0;
}
