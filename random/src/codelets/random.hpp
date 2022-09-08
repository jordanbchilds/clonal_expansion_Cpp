#ifndef __RANDOM_IPU_HPP__
#define __RANDOM_IPU_HPP__

#include <cstdint>
#include <cstddef>

#if __POPC__

#else // !__POPC__

#include <poplar/Graph.hpp>

namespace random_ipu {

  // ----------------------------------------------------------------------------
  // IPU API (as Op)
  // ----------------------------------------------------------------------------

  void random(poplar::Graph &graph,
              poplar::Tensor random,
              poplar::program::Sequence &prog);


  void setSeeds(poplar::Graph &graph,
              uint32_t base,
              uint32_t tile_increment,
              uint32_t worker_increment,
              poplar::program::Sequence &prog);

  void setRandomSeeds(poplar::Graph &graph,
              uint32_t seed,
              poplar::program::Sequence &prog);


}; // namespace random_ipu

#endif // !__POPC__

#endif // __RANDOM_IPU_HPP__
