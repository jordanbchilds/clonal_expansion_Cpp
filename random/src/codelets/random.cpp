// Common (CPU and IPU)
#include <stdint.h>
#include "random.hpp"


const size_t tiles = 1472; // TODO: Should enumerate this
const size_t workers = 6;  // TODO: Should enumerate this
const size_t per_worker = 8; // TODO: Copied from main.cpp for expediencey


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
  Vector<Input<Vector<float>>> dummy;

  bool compute(unsigned workerId) {
      const uint32_t num_workers = MultiVertex::numWorkers();
      const size_t per_worker = out.size()/num_workers;
      float* work = reinterpret_cast<float*>(&out[workerId * per_worker]);
      for (auto e = 0; e < per_worker; ++e) {
#if defined(USE_URAND_F32)
        work[e] = __builtin_ipu_urand_f32();
#elif defined(USE_URAND32)
        work[e] = __builtin_ipu_urand32();
#elif defined(USE_F32V2GRAND)
        work[e] = __builtin_ipu_f32v2grand()[1];
#else
        work[e] = dummy[workerId][e];
#endif
      }
      return true;
  }
};

#else
//!__POPC__

#include <iostream>
#include <vector>

#include <poplar/RandomSeed.hpp>

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

      // Create the dummy (pass through) tensor.
      std::vector<float>values;
      values.resize(workers * per_worker);
      for (uint32_t worker = 0; worker < workers; ++worker) {
        std::cout << "dummy worker " << worker;
        for (uint32_t e = 0; e < per_worker; ++e) {
          float v = 100 * worker + e;
          values[worker * per_worker + e] = v;
          std::cout << " " << v;
        }
        std::cout << std::endl;
      }
      auto dummy = graph.addConstant(
        poplar::FLOAT,
        {workers, per_worker},
        values.data(),
        "dummy");
      graph.connect(random["dummy"], dummy);
      for (uint32_t worker = 0; worker < workers; ++worker) {
        graph.setTileMapping(dummy, tile);
      }


      // Program sequence.
      poplar::program::Sequence seq = {
        poplar::program::Execute(random_cs)
      };

      prog.add(seq);
  }


  void setRandomSeeds(poplar::Graph &graph,
              uint32_t seed,
              poplar::program::Sequence &prog) {

    std::vector<uint32_t>values;

    // There are 4 RNG seed registers per worker per tile.
    values.resize( tiles * workers * 4 );

    std::cout << "srand " << seed << std:: endl;
    srand(seed);

    std::cout << "values size " << values.size() << std:: endl;
    for (uint32_t v = 0, tile = 0; tile < tiles; ++tile) {
      for (uint32_t worker = 0; worker < workers; ++worker) {

        uint32_t seed[4] = {
          static_cast<uint32_t>(rand()),
          static_cast<uint32_t>(rand()),
          static_cast<uint32_t>(rand()),
          static_cast<uint32_t>(rand())
        };

        if (tile == 0)
           std::cout << "v " << v << " tile " << tile << " worker " << worker << " prng " <<
             seed[0] << "," << seed[1] << "," << seed[2] << "," << seed[3] << std::endl;

        // Set the $PRNG_x_y
        // These are the xoroshiro128aox state vectors
        //  s0 = PRNG_0_1 << 32 |  PRNG_0_0
        //  s1 = PRNG_1_1 << 32 |  PRNG_0_1

        values[v++] = seed[0];
        values[v++] = seed[1];
        values[v++] = seed[2];
        values[v++] = seed[3];
      }
    }


    // Create the random seeds tensor.
    auto seeds = graph.addConstant(
        poplar::UNSIGNED_INT,
        {tiles, workers, 4},
        values.data(),
        "random_seeds");

    for (uint32_t tile = 0; tile < tiles; ++tile) {
      graph.setTileMapping(seeds[tile], tile);
    }

    // Add setSeeds op to program,
    poplar::setHwSeeds(graph, seeds, prog);
  }


  void setSeeds(poplar::Graph &graph,
              uint32_t base,
              uint32_t tile_factor,
              uint32_t worker_factor,
              poplar::program::Sequence &prog) {

    std::vector<uint32_t>values;

    // There are 4 RNG seed registers per worker per tile.
    values.resize( tiles * workers * 4 );

    std::cout << "values size " << values.size() << std:: endl;
    for (uint32_t v = 0, tile = 0; tile < tiles; ++tile) {
      for (uint32_t worker = 0; worker < workers; ++worker) {

        uint32_t seed[4] = {
          base + tile*tile_factor + worker* worker_factor,
          base + tile*tile_factor + worker* worker_factor + 1,
          base + tile*tile_factor + worker* worker_factor + 2,
          base + tile*tile_factor + worker* worker_factor + 3
        };

        if (tile == 0)
           std::cout << "v " << v << " tile " << tile << " worker " << worker << " prng " <<
             seed[0] << "," << seed[1] << "," << seed[2] << "," << seed[3] << std::endl;

        // Set the $PRNG_x_y
        // These are the xoroshiro128aox state vectors
        //  s0 = PRNG_0_1 << 32 |  PRNG_0_0
        //  s1 = PRNG_1_1 << 32 |  PRNG_0_1

        values[v++] = seed[0];
        values[v++] = seed[1];
        values[v++] = seed[2];
        values[v++] = seed[3];
      }
    }


    // Create the random seeds tensor.
    auto seeds = graph.addConstant(
        poplar::UNSIGNED_INT,
        {tiles, workers, 4},
        values.data(),
        "random_seeds");

    for (uint32_t tile = 0; tile < tiles; ++tile) {
      graph.setTileMapping(seeds[tile], tile);
    }

    // Add setSeeds op to program,
    poplar::setHwSeeds(graph, seeds, prog);
  }


} // namespace random_ipu

#endif //!__POPC__
