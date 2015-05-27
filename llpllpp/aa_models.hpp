#if ! defined(__LLPLLPLUSPLUS_AAMODEL_HPP__)
#define __LLPLLPLUSPLUS_AAMODEL_HPP__
// A simple unrooted tree. pll_utree_t in pll

#include <vector>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/model_storage_description.hpp"

namespace pllpp {

class AAModel {
  private:
  static std::vector<FixedModelRef> fixedModels;
  public:
  static const std::vector<FixedModelRef> & getFixedModels();
};

} // namespace pllpp
#endif
