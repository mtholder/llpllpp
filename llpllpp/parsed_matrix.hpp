#if ! defined(__LLPLLPLUSPLUS_PARSED_MATRIX_HPP__)
#define __LLPLLPLUSPLUS_PARSED_MATRIX_HPP__
// A simple structure storing otu name and strings for sequence data. 
// see headers and seqdata arrays in teh newick-fasta-unrooted example in PLL

#include <memory>
#include <cassert>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/otus.hpp"

namespace pllpp {

class ParsedMatrix {
  std::shared_ptr<OTUSet> otusShPtr;
  public:
  ParsedMatrix() {
  }
  static std::unique_ptr<ParsedMatrix> parseFasta(const std::string & fn,
                                                  std::shared_ptr<OTUSet> otus=nullptr);
  std::shared_ptr<OTUSet> getOTUSet() {
    return otusShPtr;
  }
};

} // namespace pllpp
#endif
