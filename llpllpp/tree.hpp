#if ! defined(__LLPLLPLUSPLUS_TREE_HPP__)
#define __LLPLLPLUSPLUS_TREE_HPP__
// A simple unrooted tree. pll_utree_t in pll

#include <memory>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/otus.hpp"

#include "pll.h"

namespace pllpp {

class UTree {
  std::shared_ptr<OTUSet> otusShPtr;
  pll_utree_t * pllTree;
  public:
  UTree()
    :pllTree(nullptr) {
  }
  ~UTree() {
    if (pllTree != nullptr) {
      pll_destroy_utree(pllTree);
    }
  }
  static std::unique_ptr<UTree> parseNewick(const std::string & fn,
                                            std::shared_ptr<OTUSet> otus=nullptr);
  void setMissingBranchLength(double edgeLength);
  std::shared_ptr<OTUSet> getOTUSet() {
    return otusShPtr;
  }
};

} // namespace pllpp
#endif
