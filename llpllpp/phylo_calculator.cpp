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
  _OperationContainer(std::size_t n,
                      std::vector<UTree::node_ptr> & traversalBuffer,
                      int traversalSize,
                      std::vector<double> &  edgeLengths,
                      std::vector<int> & matrixIndices)
    :opVec(n),
    traversalBufferRef(traversalBuffer),
    edgeLengthsRef(edgeLengths),
    matrixIndicesRef(matrixIndices) {
    createOps(traversalSize);
  }
  void createOps(int traversalSize) {
    //  given the computed traversal descriptor, generate the operations
    //    structure, and the corresponding probability matrix indices that
    //    may need recomputing 
    pll_utree_create_operations(&(traversalBufferRef[0]),
                                traversalSize,
                                &(edgeLengthsRef[0]),
                                &(matrixIndicesRef[0]),
                                ops(),
                                &matrixCount,
                                &numPendingOperations);
  }
  int getNumPendingOps() const {
    return numPendingOperations;
  }
  int getNumProbMatToCalc() const {
    return matrixCount;
  }
  pll_operation_t * ops() {
    return &(opVec[0]);
  }
  private:
  std::vector<pll_operation_t> opVec;
  std::vector<UTree::node_ptr> & traversalBufferRef;
  std::vector<double> &  edgeLengthsRef;
  std::vector<int> & matrixIndicesRef;
  int matrixCount = 0;
  int numPendingOperations = 0;
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
  tipCount = tree->getNumLeaves();
  assert(tipCount == parsedMat.getNumRows());
  probModelVec.emplace_back(msd);
  assert(tipCount > 2);
  innerNodesCount = tipCount - 2;
  nodesCount = innerNodesCount + tipCount;
  branchCount = nodesCount - 1;
  edgeLengths.resize(branchCount);
  matrixIndices.resize(branchCount);
  traversalBuffer.resize(nodesCount);
  init_traverse();
}

void PhyloCalculator::init_traverse() {
  virtualRoot = tree->pllTree;
  // allocates (if NULL) and fills:
  //    branch_lengths - lengths of postorder traversal from  node
  //    matrix_indices -- a number to each branch length, (0 to
  //       tip_count-1) for branches leading to tips, and (tip_count to
  //       2*tip_count-4) for the other branches. These numbers are used to assign
  //       slots for the probability matrix of each branch.
  //    operations - filled to compute all inner CLVs
  //assert(isInnerTernaryNode(node))
  // edge_matrix_index will point to the edge between node and node->back, 
  const auto travLen = pll_utree_traverse(virtualRoot,
                                          cb_full_traversal,
                                          &(traversalBuffer[0]));
  assert(travLen >= 0);
  opContainerPtr = new _OperationContainer(innerNodesCount,
                                           traversalBuffer,
                                           travLen,
                                           edgeLengths,
                                           matrixIndices);
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
                           opContainerPtr->getNumProbMatToCalc());

}

void PhyloCalculator::updatePartials(std::size_t ) {
  pll_update_partials(partData.partition,
                      opContainerPtr->ops(),
                      opContainerPtr->getNumPendingOps());
}

double PhyloCalculator::computeEdgeLogLikelihood(std::size_t partIndex) {
  return pll_compute_edge_loglikelihood(partData.partition, 
                                        virtualRoot->clv_index,
                                        virtualRoot->scaler_index,
                                        virtualRoot->back->clv_index,
                                        virtualRoot->back->scaler_index,
                                        virtualRoot->pmatrix_index,
                                        static_cast<int>(partIndex));

}

} // namespace
