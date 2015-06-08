#if ! defined(__LLPLLPLUSPLUS_TREE_HPP__)
#define __LLPLLPLUSPLUS_TREE_HPP__
// A simple unrooted tree. pll_utree_t in pll

#include <memory>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/otus.hpp"
struct pll_utree;
struct pll_rtree;
namespace pllpp {

template<typename W>
class WrappedTree {
  public:
  using node_type = W;
  using node_ptr = W *;
  using const_node_ptr = const W *;
  using wtree_type = WrappedTree<W>;
  private:
  std::shared_ptr<OTUSet> otusShPtr;
  node_ptr pllTree;
  WrappedTree()
    :pllTree(nullptr) {
  }
  public:
  ~WrappedTree() {
    clear();
  }
  static std::unique_ptr<WrappedTree<W> > parseNewick(const std::string & fn,
                                            std::shared_ptr<OTUSet> otus=nullptr);
  void setMissingBranchLength(double edgeLength);
  std::shared_ptr<OTUSet> getOTUSet() {
    return otusShPtr;
  }
  void clear();
  std::size_t getNumLeaves() const {
    return (otusShPtr == nullptr ? 0U : otusShPtr->size());
  }
  private:
  template<typename T> friend class PhyloCalculator;
};

using UTree = WrappedTree<pll_utree>;
using RTree = WrappedTree<pll_rtree>;


} // namespace pllpp
#endif
