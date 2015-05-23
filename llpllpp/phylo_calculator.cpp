#include "llpllpp/phylo_calculator.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {
PhyloCalculator::PhyloCalculator(const ParsedMatrix & parsedMat,
                                 const ModelStorageDescription &msd,
                                 std::shared_ptr<UTree> treePtr)
  :partData(parsedMat, msd),
  tree(treePtr),
  rateCatUpdateCounter{0},
  stateFreqUpdateCounter{0},
  exchangeUpdateCounter{0} {
  auto otus = tree->getOTUSet();
  assert(otus != nullptr);
  const int tipCount = static_cast<int>(tree->getNumLeaves());
  assert(tipCount == static_cast<int>(parsedMat.getNumRows()));
  probModelVec.emplace_back(msd);
  init_traverse();
}

void PhyloCalculator::init_traverse() {
  auto node = tree->pllTree;
  const int tipCount = static_cast<int>(tree->getNumLeaves());
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
  numPendingOperations = tipCount - 2;
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

void PhyloCalculator::updateProbMatrices(std::size_t partIndex) {
  const int partIndI = static_cast<int>(partIndex);
  const auto & model = getModel(partIndex);
  const auto sfc = model.getStateFreqCounter();
  if (stateFreqUpdateCounter.at(partIndex) != sfc) {
    const auto & sf = model.getStateFrequencies();
    pll_set_frequencies(partData.partition, partIndI, &sf[0]);
    stateFreqUpdateCounter[partIndex] = sfc;
  }
  const auto exc = model.getExchangeCounter();
  if (exchangeUpdateCounter.at(partIndex) != exc) {
    const auto & excp = model.getExchangeabilityParams();
    pll_set_subst_params(partData.partition, partIndI, &excp[0]);
    exchangeUpdateCounter[partIndex] = exc;
  }
  const auto & rateHet = model.getRateHet();
  const auto rhc = rateHet.getCounter();
  if (rateCatUpdateCounter.at(partIndex) != rhc) {
    const auto & rates = rateHet.getRates();
    pll_set_category_rates(partData.partition, &rates[0]);
    rateCatUpdateCounter[partIndex] = rhc;
  }
  pll_update_prob_matrices(partData.partition, 
                           partIndI, 
                           matrixIndices, 
                           edgeLengths, 
                           partData.numProbMats);

}

void PhyloCalculator::updatePartials(std::size_t ) {
  pll_update_partials(partData.partition, operations, numPendingOperations);
}

double PhyloCalculator::computeEdgeLogLikelihood(std::size_t partIndex) {
  return pll_compute_edge_loglikelihood(partData.partition, 
                                        clv1,
                                        clv2,
                                        edgePMatrixIndex,
                                        static_cast<int>(partIndex));

}

} // namespace
