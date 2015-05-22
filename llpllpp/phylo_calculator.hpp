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
struct pll_partition;
namespace pllpp {

class PhyloCalculator {
  PartitionedData partData;
  std::vector<DSCTProbModel> probModelVec;
  UTree & tree;
  pll_partition * partition = nullptr;
  public:
  PhyloCalculator(const ParsedMatrix & parsedMat, const ModelStorageDescription &msd, UTree & treeRef);
  DSCTProbModel & getModel(std::size_t modIndex) {
    return probModelVec.at(modIndex);
  }
  void updateProbMatrices(std::size_t partIndex);
  void updatePartials(std::size_t partIndex);
  double computeEdgeLogLikelihood(std::size_t partIndex);

};

} // namespace pllpp
#endif
