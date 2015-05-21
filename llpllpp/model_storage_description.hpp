#if ! defined(__LLPLLPLUSPLUS_MODEL_STORAGE_DESCRIPTION_HPP__)
#define __LLPLLPLUSPLUS_MODEL_STORAGE_DESCRIPTION_HPP__
// An argument to create partitions that bundles info about
//  how to allocate storage for the calculations.

#include "llpllpp/base_includes.hpp"

namespace pllpp {

class ModelStorageDescription {
  public:
  unsigned int numStates;
  unsigned int numRateCats;
  ArchAttribEnum archAttributes;
  ModelStorageDescription(unsigned int nStates, unsigned int nRateCats, ArchAttribEnum arch)
    :numStates(nStates),
    numRateCats(nRateCats),
    archAttributes(arch) {
  }
};

} // namespace pllpp
#endif
