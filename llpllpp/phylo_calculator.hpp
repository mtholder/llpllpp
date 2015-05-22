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

class PhyloCalculator {
  PartitionedData partData;
  std::vector<DSCTProbModel> probModelVec;
  UTree & tree;
  int edgePMatrixIndex;
  int clv1;
  int clv2;
  double * edgeLengths = nullptr;
  int * matrixIndices = nullptr;
  pll_operation * operations = nullptr;
  public:
  PhyloCalculator(const PhyloCalculator &) = delete;
  PhyloCalculator & operator=(const PhyloCalculator &) = delete;
  PhyloCalculator(const ParsedMatrix & parsedMat, const ModelStorageDescription &msd, UTree & treeRef);
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
