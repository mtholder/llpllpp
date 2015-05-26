#if ! defined(__LLPLLPLUSPLUS_TREE_HPP__)
#define __LLPLLPLUSPLUS_TREE_HPP__
// A simple unrooted tree. pll_utree_t in pll

#include <memory>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/otus.hpp"
struct pll_utree;
namespace pllpp {

class UTree {
  std::shared_ptr<OTUSet> otusShPtr;
  pll_utree * pllTree;
  UTree()
    :pllTree(nullptr) {
  }
  public:
  ~UTree() {
    clear();
  }
  static std::unique_ptr<UTree> parseNewick(const std::string & fn,
                                            std::shared_ptr<OTUSet> otus=nullptr);
  void setMissingBranchLength(double edgeLength);
  std::shared_ptr<OTUSet> getOTUSet() {
    return otusShPtr;
  }
  void clear();
  std::size_t getNumLeaves() const {
    return (otusShPtr == nullptr ? 0U : otusShPtr->size());
  }
  friend class PhyloCalculator;
};

} // namespace pllpp
#endif
