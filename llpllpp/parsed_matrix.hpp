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
  std::vector<const char *> headers;
  std::vector<const char *> seqData;
  std::size_t seqLen = 0U;
  public:
  ParsedMatrix(std::shared_ptr<OTUSet> otus)
    :otusShPtr(otus) {
  }
  ~ParsedMatrix() {
    clear();
  }
  static std::unique_ptr<ParsedMatrix> parseFasta(const std::string & fn,
                                                  std::shared_ptr<OTUSet> otus=nullptr);
  std::shared_ptr<OTUSet> getOTUSet() const {
    return otusShPtr;
  }
  void clear();
  std::size_t getLength() const {
    return seqLen;
  }
  std::size_t getNumRows() const {
    assert(headers.size() == seqData.size());
    return headers.size();
  }
  const char * getName(std::size_t index) const {
    return headers.at(index);
  }
  const char * getData(std::size_t index) const {
    return seqData.at(index);
  }
};

} // namespace pllpp
#endif
