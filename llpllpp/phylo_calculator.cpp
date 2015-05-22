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
  init_traverse();
}

void PhyloCalculator::init_traverse() {
  auto node = tree.pllTree;
  const int tipCount = static_cast<int>(tree.getNumLeaves());
  // allocates (if NULL) and fills:
  //    branch_lengths - lengths of postorder traversal from  node
  //    matrix_indices -- a number to each branch length, (0 to
  //       tip_count-1) for branches leading to tips, and (tip_count to
  //       2*tip_count-4) for the other branches. These numbers are used to assign
  //       slots for the probability matrix of each branch.
  //    operations - filled to compute all inner CLVs
  //assert(isInnerTernaryNode(node))
  // edge_matrix_index will point to the edge between node and node->back, 
  // clv1 is set to the index of node->back and 
  // clv2 to the index of node. 
  // edge_matrix_index, clv1, and clv2 can be used to evaluate the log-likelihood using the
  //   pll_compute_edge_loglikelihood function 
  pll_traverse_utree(node,
                     tipCount,
                     &edgeLengths,
                     &matrixIndices,
                     &operations,
                     &edgePMatrixIndex,
                     &clv1, 
                     &clv2);
}

void PhyloCalculator::clear() {
  partData.clear();
}

} // namespace
