#include "llpllpp/tree.hpp"
#include "pll.h"
namespace pllpp {

template<typename T>
typename T::node_ptr _rawRead(const std::string & rn, int & tipCount);
template<typename T>
std::size_t _queryTipnodes(typename T::node_ptr, typename T::node_ptr *);
template<typename T>
std::size_t _queryInternals(typename T::node_ptr, typename T::node_ptr *);
template<typename T>
void _destroyTree(typename T::node_ptr & tree);
template<typename T>
void set_missing_branch_length_recursive(T * tree, double length);
template<typename T>
T * _firstRecurseTop(T * tree);
template<typename T>
T * _secondRecurseTop(T * tree);

template<>
inline UTree::node_ptr _rawRead<UTree>(const std::string & fn,
                                int & tipCount) {
  return pll_utree_parse_newick(fn.c_str(), &tipCount);
}

template<>
inline RTree::node_ptr _rawRead<RTree>(const std::string & fn,
                                int & tipCount) {
  return pll_rtree_parse_newick(fn.c_str(), &tipCount);
}

template<>
inline std::size_t _queryTipnodes<UTree>(UTree::node_ptr t,
                                 UTree::node_ptr * v) {
  return static_cast<std::size_t>(pll_utree_query_tipnodes(t, v));
}

template<>
inline std::size_t _queryInternals<UTree>(UTree::node_ptr t,
                                  UTree::node_ptr * v) {
  return static_cast<std::size_t>(pll_utree_query_innernodes(t, v));
}

template<>
inline std::size_t _queryTipnodes<RTree>(RTree::node_ptr t,
                          RTree::node_ptr * v) {
  return static_cast<std::size_t>(pll_rtree_query_tipnodes(t, v));
}

template<>
inline std::size_t _queryInternals<RTree>(RTree::node_ptr t,
                                  RTree::node_ptr * v) {
  return static_cast<std::size_t>(pll_rtree_query_innernodes(t, v));
}


template<>
inline void _destroyTree<UTree>(UTree::node_ptr & tree) {
  pll_utree_destroy(tree);
  tree = nullptr;
}

template<>
inline void _destroyTree<RTree>(RTree::node_ptr & tree) {
  pll_rtree_destroy(tree);
  tree = nullptr;
}


template<>
inline pll_utree * _firstRecurseTop<pll_utree>(pll_utree * tree) {
  return tree;
}

template<>
inline pll_utree * _secondRecurseTop<pll_utree>(pll_utree * tree) {
  return tree->back;
}

template<>
inline pll_rtree * _firstRecurseTop<pll_rtree>(pll_rtree * tree) {
  return tree->left;
}

template<>
inline pll_rtree * _secondRecurseTop<pll_rtree>(pll_rtree * tree) {
  return tree->right;
}


// from nfu.c
template<>
void set_missing_branch_length_recursive<pll_utree>(pll_utree * tree, 
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
    set_missing_branch_length_recursive<pll_utree>(tree->next->back, length);
    set_missing_branch_length_recursive<pll_utree>(tree->next->next->back, length);
  }
}

// from nfu.c
template<>
void set_missing_branch_length_recursive<pll_rtree>(pll_rtree * tree, 
                                                    double length) {
  if (tree == nullptr) {
    return;
  }
  if (tree->length <= 0.0) {
    tree->length = length;
  }

  if (tree->left)
    set_missing_branch_length_recursive<pll_rtree>(tree->left, length);

  if (tree->right)
    set_missing_branch_length_recursive<pll_rtree>(tree->right, length);
}


template<typename W>
std::unique_ptr<WrappedTree<W> >
WrappedTree<W>::parseNewick(const std::string & fn,
                                          std::shared_ptr<OTUSet> otus) {
  int tipCount;
  node_ptr tree = _rawRead<wtree_type>(fn, tipCount);
  const std::size_t tipCountU = static_cast<std::size_t>(tipCount);
  if (tree == nullptr) {
    throw PLLException(std::string("Could not read a tree from ") + fn);
  }
  wtree_type * utree = new wtree_type();
  utree->nodes.clear();
  utree->nodes.resize(2*tipCountU);
  try {
    if (otus == nullptr) {
      otus = std::make_shared<OTUSet>();
      const auto qtc = _queryTipnodes<wtree_type>(tree, &(utree->nodes[0]));
      assert(qtc == tipCountU);
      for (auto i = 0UL ; i < tipCountU; ++i) {
        assert(utree->nodes[i] != nullptr);
        auto j = otus->addNewName(utree->nodes[i]->label);
        assert(j == static_cast<std::size_t>(i));
      }
    } else {
      NOT_IMPLEMENTED
    }
  } catch (...) {
    _destroyTree<wtree_type>(tree);
    delete utree;
    throw;
  }
  const auto nInternals = _queryInternals<wtree_type>(tree, &(utree->nodes[tipCountU]));
  assert(nInternals + tipCountU <= utree->nodes.size());
  utree->nodes.resize(nInternals + tipCountU);
  utree->pllTree = tree;
  utree->otusShPtr = otus;
  utree->_initEdges();
  return std::unique_ptr<wtree_type>(utree);
}

// fills `edges` based on other fields. callec by parseNewick
template<>
void UTree::_initEdges() {
  assert(edges.empty());
  assert(!nodes.empty());
}

// fills `edges` based on other fields. callec by parseNewick
template<>
void RTree::_initEdges() {
  assert(edges.empty());
  assert(!nodes.empty());
  edges.resize(nodes.size() - 1);
}

template<typename W>
void WrappedTree<W>::setMissingBranchLength(double edgeLength) {
  assert(pllTree != nullptr);
  set_missing_branch_length_recursive<node_type>(_firstRecurseTop<node_type>(pllTree), edgeLength);
  set_missing_branch_length_recursive<node_type>(_secondRecurseTop<node_type>(pllTree), edgeLength);
}

template<typename W>
void WrappedTree<W>::clear() {
  if (pllTree != nullptr) {
   _destroyTree<wtree_type>(pllTree);
  }
}

template class WrappedTree<pll_utree>; // force explicit instantiaion of this template.
template class WrappedTree<pll_rtree>; // force explicit instantiaion of this template.

} // namespace
