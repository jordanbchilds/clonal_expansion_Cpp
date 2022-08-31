// Common (CPU and IPU)
#include <stdint.h>
#include "random.hpp"

#ifdef __POPC__

// IPU-specific code.

#include <poplar/Vertex.hpp>

using namespace poplar;

// Codelets to decode a JPEG file.

// Random
class Random : public MultiVertex {
public:
  Random();
  Output<Vector<float>> out;

  bool compute(unsigned workerId) {
      const uint32_t num_workers = MultiVertex::numWorkers();
      const size_t per_worker = out.size()/num_workers;
      float* work = reinterpret_cast<float*>(&out[workerId * per_worker]);
	  
      for (auto e = 0; e < per_worker; ++e) {
			  work[e] = __builtin_ipu_urand_f32();
      }
      return true;
  }
};

#else
//!__POPC__

#include <iostream>
#include <vector>

namespace random_ipu {

  // ----------------------------------------------------------------------------
  // IPU API (as Op)
  // ----------------------------------------------------------------------------
  void random(poplar::Graph &graph,
              poplar::Tensor random_out,
              poplar::program::Sequence &prog) {

      // TODO:
      // Select best tile from tensor location
      const uint32_t tile = 0;

      graph.setTileMapping(random_out, tile);

      graph.addCodelets("random.gp");

      // Each compute set.
      auto random_cs = graph.addComputeSet("Random");

      // Add vertices.
      auto random = graph.addVertex(random_cs, "Random");

      // Map vertices to our tile and connect common buffers.
      graph.setTileMapping(random, tile);
      graph.connect(random["out"], random_out);

      // Program sequence.
      poplar::program::Sequence seq = {
        poplar::program::Execute(random_cs)
      };
      prog.add(seq);
  }

} // namespace random_ipu

#endif //!__POPC__
