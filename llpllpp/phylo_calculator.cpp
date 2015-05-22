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
  assert(otus != nullptr);
  const int tipCount = static_cast<int>(tree.getNumLeaves());
  assert(tipCount == static_cast<int>(parsedMat.getNumRows()));
  probModelVec.emplace_back(msd);
}

void PhyloCalculator::clear() {
  partData.clear();
}

} // namespace
