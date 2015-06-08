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

template<typename T> class _OperationContainer;

template<typename W>
class PhyloCalculator {
  using node_type = typename W::node_type;
  using node_ptr = typename W::node_ptr;
  using const_node_ptr = const node_type *;
  PartitionedData partData;
  std::vector<DSCTProbModel> probModelVec;
  std::shared_ptr<W> tree;
  const_node_ptr virtualRoot = nullptr;
  std::vector<double> edgeLengths;
  std::vector<int> matrixIndices;
  _OperationContainer<node_type> * opContainerPtr;
  std::vector<node_ptr> traversalBuffer;
  typedef std::vector<unsigned long> UpdateCounterVec;
  UpdateCounterVec rateCatUpdateCounter;
  UpdateCounterVec stateFreqUpdateCounter;
  UpdateCounterVec exchangeUpdateCounter;
  // convenience. copies of limits..
  std::size_t tipCount = 0U;
  std::size_t innerNodesCount = 0U;
  std::size_t nodesCount = 0U;
  std::size_t branchCount = 0U;
  mutable std::vector<node_ptr> innerNodesAliases;
  public:
  PhyloCalculator(const PhyloCalculator &) = delete;
  PhyloCalculator & operator=(const PhyloCalculator &) = delete;
  PhyloCalculator(const ParsedMatrix & parsedMat,
                  const ModelStorageDescription &msd,
                  std::shared_ptr<W> treeRef);
  ~PhyloCalculator() {
    clear();
  }
  void clear();
  DSCTProbModel & getModel(std::size_t modIndex) {
    return probModelVec.at(modIndex);
  }
  void updateProbMatrices(std::size_t partIndex);
  void updatePartials(std::size_t partIndex);
  double computeLogLikelihood(std::size_t partIndex);
  std::size_t getInnerNodesCount() const {
    return innerNodesCount;
  }
  const W & getTree() const {
    assert(tree != nullptr);
    return *tree;
  }
  void setVirtualRoot(const_node_ptr v) {
    partialTraverse(v);
  }
  node_ptr * getInnerNodesArray() const {
    if (innerNodesAliases.empty()) {
      innerNodesAliases.resize(getInnerNodesCount());
    }
    _fillInnerNodesArray(&(innerNodesAliases[0]));
    return &(innerNodesAliases[0]);
  }
  private:
  void partialTraverse(const_node_ptr v);
  void initTraverse();
  void _fillInnerNodesArray(node_ptr * arr) const;
};

using UPhyloCalculator = PhyloCalculator<UTree>;
using RPhyloCalculator = PhyloCalculator<RTree>;

} // namespace pllpp
#endif
