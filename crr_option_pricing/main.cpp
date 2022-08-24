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
  bool isAmericanBool = true;
  bool isCallBool = true;

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
  graph.addCodelets("compute_codelet.cpp");

  // Create a control program that is a sequence of steps
  Sequence prog;

  int nVals[datasetSize];
  float SVals[datasetSize];
  float KVals[datasetSize];
  float rVals[datasetSize];
  float qVals[datasetSize];
  float vVals[datasetSize];
  float tVals[datasetSize];
  bool callVals[datasetSize];
  bool americanVals[datasetSize];
  for (int i = 0; i < datasetSize; ++i){
    nVals[i] = nVal;
    SVals[i] = 100 + 0.01 * i;
    KVals[i] = 100;
    rVals[i] = 0.01;
    qVals[i] = 0;
    vVals[i] = 0.5;
    tVals[i] = 0.1;
    callVals[i] = isCallBool;
    americanVals[i] = isAmericanBool;
  }

  // Add steps to initialize the variables
  Tensor n = graph.addConstant<int>(INT, {datasetSize}, nVals);
  Tensor S = graph.addConstant<float>(FLOAT, {datasetSize}, SVals);
  Tensor K = graph.addConstant<float>(FLOAT, {datasetSize}, KVals);
  Tensor r = graph.addConstant<float>(FLOAT, {datasetSize}, rVals);
  Tensor q = graph.addConstant<float>(FLOAT, {datasetSize}, qVals);
  Tensor v = graph.addConstant<float>(FLOAT, {datasetSize}, vVals);
  Tensor T = graph.addConstant<float>(FLOAT, {datasetSize}, tVals);
  Tensor isCallOption = graph.addConstant<bool>(BOOL, {datasetSize}, callVals);
  Tensor isAmerican = graph.addConstant<bool>(BOOL, {datasetSize}, americanVals);
  long unsigned int vecSize=nVal + 1;
  Tensor prices = graph.addVariable(FLOAT, {datasetSize, vecSize}, "prices");
  Tensor output = graph.addVariable(FLOAT, {datasetSize}, "output");

  ComputeSet computeSet = graph.addComputeSet("computeSet");

  // iterate through tiles on the IPU, map 6 sets of variables and 6 option pricers to each tile
  for (int i = 0; i < datasetSize; ++i)
  {
    int roundCount = i % int(numberOfTiles * threadsPerTile);
    int tileInt = std::floor(
        float(roundCount) / float(threadsPerTile)
        );
    // std::cout << "theadsPerIPU=" << int(numberOfTiles * threadsPerTile) << "i="<< i << " roundCount=" << roundCount << " tileInt=" << tileInt << std::endl;

    graph.setTileMapping(n[i], tileInt);
    graph.setTileMapping(S[i], tileInt);
    graph.setTileMapping(K[i], tileInt);
    graph.setTileMapping(r[i], tileInt);
    graph.setTileMapping(q[i], tileInt);
    graph.setTileMapping(v[i], tileInt);
    graph.setTileMapping(T[i], tileInt);
    graph.setTileMapping(isCallOption[i], tileInt);
    graph.setTileMapping(isAmerican[i], tileInt);
    graph.setTileMapping(prices[i], tileInt);
    graph.setTileMapping(output[i], tileInt);

    VertexRef vtx = graph.addVertex(computeSet, "CRRVertex");
    graph.setTileMapping(vtx, tileInt);

    graph.connect(vtx["n"], n[i]);
    graph.connect(vtx["S"], S[i]);
    graph.connect(vtx["K"], K[i]);
    graph.connect(vtx["r"], r[i]);
    graph.connect(vtx["q"], q[i]);
    graph.connect(vtx["v"], v[i]);
    graph.connect(vtx["T"], T[i]);
    graph.connect(vtx["isCallOption"], isCallOption[i]);
    graph.connect(vtx["isAmerican"], isAmerican[i]);
    graph.connect(vtx["prices"], prices[i]);
    graph.connect(vtx["out"], output[i]);

  }

  // Add a step to execute the compute set
  prog.add(Execute(computeSet));
  // Add a step to print out the option values
  prog.add(PrintTensor("output", output));
  // Create the engine
  Engine engine(graph, prog);
  engine.load(device);

  auto start = std::chrono::system_clock::now();
  // Run the control program
  engine.run(0);
  auto end = std::chrono::system_clock::now();

  std::cout << "Processed " << datasetSize << " options" << std::endl;
  std::cout << "N=" << nVal << std::endl;
  std::cout << "isAmerican=" << isAmericanBool << std::endl;

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::cout << "Completed computation at " << std::ctime(&end_time)
            << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
  
  std::cout << "Rate=" << (float(datasetSize)/elapsed_seconds.count()) << std::endl;

  return 0;
}
