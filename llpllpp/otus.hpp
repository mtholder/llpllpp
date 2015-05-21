#if ! defined(__LLPLLPLUSPLUS_OTUS_HPP__)
#define __LLPLLPLUSPLUS_OTUS_HPP__
/// holds mapping of tip label to integer
#include <map>
#include "llpllpp/base_includes.hpp"

namespace pllpp {

class OTUSet {
  std::map<std::string, std::size_t> nameToIndex;
  std::vector<const char *> names;
  public:
  std::size_t size() const {
    return nameToIndex.size();
  }
  // throws PLLException if n has already been stored...
  std::size_t addNewName(const std::string &n) {
    std::size_t ind = names.size();
    auto itAddedPair = nameToIndex.emplace(n, ind);
    if (!itAddedPair.second) {
      std::string message = std::string("Illegal of tip name \"") + n + std::string("\"");
      throw PLLException(message);
    }
    names.push_back(itAddedPair.first->first.c_str());
    return ind;
  }
};

} // namespace pllpp
#endif
