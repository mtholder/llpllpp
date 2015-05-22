#include "llpllpp/phylo_calculator.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {
PhyloCalculator::PhyloCalculator(const ParsedMatrix & parsedMat,
                                 const ModelStorageDescription &msd,
                                 UTree & treeRef)
  :partData(parsedMat, msd),
  tree(treeRef) {
  auto otus = tree.getOTUSet();
  if (otus == nullptr) {
    throw PLLException("Programmer Error: tree must have otuSet before calling PhyloCalculator ctor\n");
  }
  const int tipCount = static_cast<int>(tree.getNumLeaves());
  probModelVec.emplace_back(msd);
  const auto numModels = 1; // for this ctor only
  const auto numProbMats = 2*tipCount - 3;
  partition = pll_create_partition(tipCount,
                                   tipCount - 2,
                                   static_cast<int>(msd.numStates),
                                   static_cast<int>(parsedMat.getLength()),
                                   numModels,
                                   numProbMats,
                                   static_cast<int>(msd.numRateCats),
                                   tipCount - 2,
                                   PLL_ATTRIB_ARCH_SSE);
}

} // namespace
