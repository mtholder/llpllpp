#include "llpllpp/tree.hpp"
#include "pll.h"

namespace pllpp {
std::unique_ptr<UTree> UTree::parseNewick(const std::string & fn,
                                          std::shared_ptr<OTUSet> otus) {
  int tipCount;
  pll_utree_t * tree = pll_utree_parse_newick(fn.c_str(), &tipCount);
  if (tree == nullptr) {
    throw PLLException(std::string("Could not read a tree from ") + fn);
  }
  // Move
  std::vector<node_ptr> tipNodes;
  tipNodes.resize(static_cast<std::size_t>(tipCount));

  try {
    if (otus == nullptr) {
      otus = std::make_shared<OTUSet>();
      pll_utree_query_tipnodes(tree, &(tipNodes[0]));
      for (auto i = 0UL ; i < static_cast<std::size_t>(tipCount); ++i) {
        assert(tipNodes[i] != nullptr);
        auto j = otus->addNewName(tipNodes[i]->label);
        assert(j == static_cast<std::size_t>(i));
      }
    } else {
      NOT_IMPLEMENTED
    }
  } catch (...) {
    pll_utree_destroy(tree);
    throw;
  }
  // Now transfer to a new UTree

  UTree * utree = new UTree();
  utree->pllTree = tree;
  utree->otusShPtr = otus;
  return std::unique_ptr<UTree>(utree);
}

// from nfu.c
static void set_missing_branch_length_recursive(pll_utree_t * tree, 
                                                double length) {
  if (tree == nullptr) {
    return;
  }
  /* set branch length to default if not set */
  if (tree->length <= 0.0) {
      tree->length = length;
  }
  if (tree->next) {
    if (tree->next->length <= 0.0) {
      tree->next->length = length;
    }
    if (tree->next->next->length <= 0.0) {
      tree->next->next->length = length;
    }
    set_missing_branch_length_recursive(tree->next->back, length);
    set_missing_branch_length_recursive(tree->next->next->back, length);
  }
}

void UTree::setMissingBranchLength(double edgeLength) {
  set_missing_branch_length_recursive(pllTree, edgeLength);
  set_missing_branch_length_recursive(pllTree->back, edgeLength);
}

void UTree::clear() {
  if (pllTree != nullptr) {
    pll_utree_destroy(pllTree);
    pllTree = nullptr;
  }
}

} // namespace
