#if ! defined(__LLPLLPLUSPLUS_PARTITIONED_DATA_HPP__)
#define __LLPLLPLUSPLUS_PARTITIONED_DATA_HPP__
// The data and calc buffer storage for calculations. Similar to pll_partition_t
#include <memory>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/otus.hpp"
#include "llpllpp/parsed_matrix.hpp"
#include "llpllpp/model_storage_description.hpp"
#include "llpllpp/dsct_prob_model.hpp"
struct pll_partition;
namespace pllpp {

class PartitionedData {
  pll_partition * partition = nullptr;
  int numProbMats = 0;
  public:
  PartitionedData(const ParsedMatrix & parsedMat, const ModelStorageDescription &);
  ~PartitionedData() {
    clear();
  }
  void clear();
  template<typename T> friend class PhyloCalculator;
};

} // namespace pllpp
#endif
