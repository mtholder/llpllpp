#include "llpllpp/phylo_calculator.hpp"
#include "pll.h"
#include <cstdlib>
extern "C" {
/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_utree_t * node);

static int cb_full_traversal(pll_utree_t *)
{
  return 1;
}

} // end extern C

namespace pllpp {

class _OperationContainer {
  public:
  _OperationContainer(std::size_t n)
    :opVec(n) {
  }
  pll_operation_t * ops() {
    return &(opVec[0]);
  }
  private:
  std::vector<pll_operation_t> opVec;
  friend class PhyloCalculator;
};

PhyloCalculator::PhyloCalculator(const ParsedMatrix & parsedMat,
                                 const ModelStorageDescription &msd,
                                 std::shared_ptr<UTree> treePtr)
  :partData(parsedMat, msd),
  tree(treePtr),
  opContainerPtr(nullptr),
  rateCatUpdateCounter{0},
  stateFreqUpdateCounter{0},
  exchangeUpdateCounter{0} {
  auto otus = tree->getOTUSet();
  assert(otus != nullptr);
  const int tipCount = static_cast<int>(tree->getNumLeaves());
  assert(tipCount == static_cast<int>(parsedMat.getNumRows()));
  probModelVec.emplace_back(msd);
  assert(tipCount > 2);
  const std::size_t innerNodesCount = static_cast<std::size_t>(tipCount) - 2;
  const std::size_t nodesCount = innerNodesCount + static_cast<std::size_t>(tipCount);
  const std::size_t branchCount = nodesCount - 1;
  edgeLengths.resize(branchCount);
  matrixIndices.resize(branchCount);
  opContainerPtr = new _OperationContainer(innerNodesCount);
  traversalBuffer.resize(nodesCount);
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
  const auto rc = pll_utree_traverse(node,
                                    cb_full_traversal,
                                    &(traversalBuffer[0]));
  assert(rc != -1);
}

void PhyloCalculator::clear() {
  partData.clear();
  if (opContainerPtr) {
    delete opContainerPtr;
    opContainerPtr = nullptr;
  }
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
                           &(matrixIndices[0]), 
                           &(edgeLengths[0]), 
                           partData.numProbMats);

}

void PhyloCalculator::updatePartials(std::size_t ) {
  pll_update_partials(partData.partition, opContainerPtr->ops(), numPendingOperations);
}

double PhyloCalculator::computeEdgeLogLikelihood(std::size_t partIndex) {
  return pll_compute_edge_loglikelihood(partData.partition, 
                                        clv1,
                                        scaler1Index,
                                        clv2,
                                        scaler2Index,
                                        edgePMatrixIndex,
                                        static_cast<int>(partIndex));

}

} // namespace
