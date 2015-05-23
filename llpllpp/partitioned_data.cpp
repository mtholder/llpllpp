#include "llpllpp/partitioned_data.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {
PartitionedData::PartitionedData(const ParsedMatrix & parsedMat,
                                 const ModelStorageDescription &msd) {
  auto otus = parsedMat.getOTUSet();
  assert(otus != nullptr);
  const int tipCount = static_cast<int>(otus->size());
  assert(tipCount == static_cast<int>(parsedMat.getNumRows()));
  const auto numModels = 1; // for this ctor only
  numProbMats = 2*tipCount - 3;
  const auto numScaleBuffers = tipCount - 2;
  partition = pll_create_partition(tipCount,
                                   tipCount - 2,
                                   static_cast<int>(msd.numStates),
                                   static_cast<int>(parsedMat.getLength()),
                                   numModels,
                                   numProbMats,
                                   static_cast<int>(msd.numRateCats),
                                   numScaleBuffers,
                                   PLL_ATTRIB_ARCH_SSE);
  if (partition == nullptr) {
    throw PLLException("Could not allocate a PartitionedData partition data structure");
  }
  for (std::size_t i = 0; i < static_cast<std::size_t>(tipCount); ++i) {
    const auto name = parsedMat.getName(i);
    const int tip_clv_index = static_cast<int>(otus->getIndex(name));
    if (pll_set_tip_states(partition, tip_clv_index, pll_map_nt, parsedMat.getData(i)) == PLL_FAILURE) {
      throw PLLException(std::string("Problem setting tip for \"") + name + std::string("\""));
    }
  }

}

void PartitionedData::clear() {
  if (partition != nullptr) {
    pll_destroy_partition(partition);
    partition = nullptr;
  }
}

} // namespace
