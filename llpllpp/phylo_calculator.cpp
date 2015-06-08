#include "llpllpp/phylo_calculator.hpp"
#include "pll.h"
#include <cstdlib>
extern "C" {
/* a callback function for performing a full traversal */
static int cb_full_utraversal(pll_utree_t * node);
static int cb_full_rtraversal(pll_rtree_t * node);

static int cb_full_rtraversal(pll_rtree_t *)
{
  return 1;
}

static int cb_full_utraversal(pll_utree_t *)
{
  return 1;
}

} // end extern C

namespace pllpp {
template<typename T>
class _OperationContainer {
  public:
  _OperationContainer(std::size_t n,
                      std::vector<T *> & traversalBuffer,
                      int traversalSize,
                      std::vector<double> &  edgeLengths,
                      std::vector<int> & matrixIndices)
    :opVec(n),
    traversalBufferRef(traversalBuffer),
    edgeLengthsRef(edgeLengths),
    matrixIndicesRef(matrixIndices) {
    createOps(traversalSize);
  }
  void createOps(int traversalSize);
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
  std::vector<T *> & traversalBufferRef;
  std::vector<double> &  edgeLengthsRef;
  std::vector<int> & matrixIndicesRef;
  int matrixCount = 0;
  int numPendingOperations = 0;
  template<typename W> friend class PhyloCalculator;
};

template <>
inline void _OperationContainer<pll_utree_t>::createOps(int traversalSize) {
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

template <>
inline void _OperationContainer<pll_rtree_t>::createOps(int traversalSize) {
    //  given the computed traversal descriptor, generate the operations
    //    structure, and the corresponding probability matrix indices that
    //    may need recomputing 
    pll_rtree_create_operations(&(traversalBufferRef[0]),
                                traversalSize,
                                &(edgeLengthsRef[0]),
                                &(matrixIndicesRef[0]),
                                ops(),
                                &matrixCount,
                                &numPendingOperations);
}

template<typename T>
int _callPLLTraverse(T * node, int cb(T*), T ** travBuff);

template<typename T>
int _callFullPLLTraverse(T * node, T ** travBuff);

template<>
inline int _callPLLTraverse<pll_utree_t>(pll_utree_t * node, int cb(pll_utree_t *), pll_utree_t ** travBuff) {
  return pll_utree_traverse(node, cb, travBuff);
}

template<>
inline int _callFullPLLTraverse<pll_utree_t>(pll_utree_t * node, pll_utree_t ** travBuff) {
  return _callPLLTraverse<pll_utree_t>(node,  cb_full_utraversal, travBuff);
}

template<>
inline int _callPLLTraverse<pll_rtree_t>(pll_rtree_t * node, int cb(pll_rtree_t *), pll_rtree_t ** travBuff) {
  return pll_rtree_traverse(node, cb, travBuff);
}

template<>
inline int _callFullPLLTraverse<pll_rtree_t>(pll_rtree_t * node, pll_rtree_t ** travBuff) {
  return _callPLLTraverse<pll_rtree_t>(node,  cb_full_rtraversal, travBuff);
}


template<typename W>
PhyloCalculator<W>::PhyloCalculator(const ParsedMatrix & parsedMat,
                                    const ModelStorageDescription &msd,
                                    std::shared_ptr<W> treePtr)
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

template<typename W>
void PhyloCalculator<W>::init_traverse() {
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
  const auto travLen = _callFullPLLTraverse<node_type>(virtualRoot,
                                                       &(traversalBuffer[0]));
  assert(travLen >= 0);
  opContainerPtr = new _OperationContainer<node_type>(innerNodesCount,
                                                      traversalBuffer,
                                                      travLen,
                                                      edgeLengths,
                                                      matrixIndices);
}

template<typename W>
void PhyloCalculator<W>::clear() {
  partData.clear();
  if (opContainerPtr) {
    delete opContainerPtr;
    opContainerPtr = nullptr;
  }
}

template<typename W>
void PhyloCalculator<W>::updateProbMatrices(std::size_t partIndex) {
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

template<typename W>
void PhyloCalculator<W>::updatePartials(std::size_t ) {
  pll_update_partials(partData.partition,
                      opContainerPtr->ops(),
                      opContainerPtr->getNumPendingOps());
}

template<>
double PhyloCalculator<UTree>::computeLogLikelihood(std::size_t partIndex) {
  return pll_compute_edge_loglikelihood(partData.partition, 
                                        virtualRoot->clv_index,
                                        virtualRoot->scaler_index,
                                        virtualRoot->back->clv_index,
                                        virtualRoot->back->scaler_index,
                                        virtualRoot->pmatrix_index,
                                        static_cast<int>(partIndex));
}

template<>
double PhyloCalculator<RTree>::computeLogLikelihood(std::size_t partIndex) {
  return pll_compute_root_loglikelihood(partData.partition, 
                                        virtualRoot->clv_index,
                                        virtualRoot->scaler_index,
                                        static_cast<int>(partIndex));
}
template class PhyloCalculator<UTree>;
template class PhyloCalculator<RTree>;

} // namespace
