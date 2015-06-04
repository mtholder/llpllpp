#if ! defined(__LLPLLPLUSPLUS_PHYLO_CALCULATOR_HPP__)
#define __LLPLLPLUSPLUS_PHYLO_CALCULATOR_HPP__
// Acts like wrapper around tree and PartitionData. similar to the partition in PLL

#include <memory>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/model_storage_description.hpp"
#include "llpllpp/parsed_matrix.hpp"
#include "llpllpp/partitioned_data.hpp"
#include "llpllpp/tree.hpp"
struct pll_operation;
namespace pllpp {

class _OperationContainer;

class PhyloCalculator {
  PartitionedData partData;
  std::vector<DSCTProbModel> probModelVec;
  std::shared_ptr<UTree> tree;
  UTree::node_ptr virtualRoot = nullptr;
  std::vector<double> edgeLengths;
  std::vector<int> matrixIndices;
  _OperationContainer * opContainerPtr;
  std::vector<UTree::node_ptr> traversalBuffer;
  typedef std::vector<unsigned long> UpdateCounterVec;
  UpdateCounterVec rateCatUpdateCounter;
  UpdateCounterVec stateFreqUpdateCounter;
  UpdateCounterVec exchangeUpdateCounter;
  // convenience. copies of limits..
  std::size_t tipCount = 0U;
  std::size_t innerNodesCount = 0U;
  std::size_t nodesCount = 0U;
  std::size_t branchCount = 0U;
  public:
  PhyloCalculator(const PhyloCalculator &) = delete;
  PhyloCalculator & operator=(const PhyloCalculator &) = delete;
  PhyloCalculator(const ParsedMatrix & parsedMat,
                  const ModelStorageDescription &msd,
                  std::shared_ptr<UTree> treeRef);
  ~PhyloCalculator() {
    clear();
  }
  void clear();
  DSCTProbModel & getModel(std::size_t modIndex) {
    return probModelVec.at(modIndex);
  }
  void updateProbMatrices(std::size_t partIndex);
  void updatePartials(std::size_t partIndex);
  double computeEdgeLogLikelihood(std::size_t partIndex);
  private:
  void init_traverse();

};

} // namespace pllpp
#endif
