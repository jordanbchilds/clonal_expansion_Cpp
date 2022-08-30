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


}; // namespace random_ipu

#endif // !__POPC__

#endif // __RANDOM_IPU_HPP__
