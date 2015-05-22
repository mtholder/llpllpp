#include "llpllpp/parsed_matrix.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {

std::unique_ptr<ParsedMatrix> ParsedMatrix::parseFasta(const std::string & fn,
                                                       std::shared_ptr<OTUSet> otus) {
  pll_fasta_t * fp = pll_fasta_open(fn.c_str(), pll_map_fasta);
  if (!fp){
    throw PLLException(std::string("Error opening file \"") + fn + std::string("\""));
  }
  assert(otus != nullptr);
  const auto tipCount = otus->size();
  std::unique_ptr<ParsedMatrix> pm = std::make_unique<ParsedMatrix>(otus);
  pm->headers.resize(tipCount);
  pm->seqData.resize(tipCount);
  long hdrlen;
  long seqno;
  long currSeqLen;
  auto i = 0U;
  char * hdr = nullptr;
  char * seq = nullptr;
  try {
    for (; pll_fasta_getnext(fp, &hdr, &hdrlen, &seq, &currSeqLen, &seqno); ++i) {
      if (i >= tipCount) {
        throw PLLException("FASTA file contains more sequences than expected");
      }
      if ((pm->seqLen != 0 && pm->seqLen != static_cast<std::size_t>(currSeqLen))
          || currSeqLen < 0) {
        throw PLLException(std::string("FASTA file does not contain equal size sequences") 
                           + std::to_string(i + 1) + std::string(" with header \"")
                           + std::string(hdr) + std::string("\".\n"));
      }
      pm->seqLen = static_cast<std::size_t>(currSeqLen);
      pm->headers[i] = hdr;
      pm->seqData[i] = seq;
    }
      /* did we stop reading the file because we reached EOF? */
    if (pll_errno != PLL_ERROR_FILE_EOF) {
      throw PLLException(std::string("Error while reading FASTA file \"")
                        + fn + std::string("\""));
    }
  } catch(...) {
    pll_fasta_close(fp);
    throw;
  }
  pll_fasta_close(fp);
  if (pm->seqLen == 0){
      throw PLLException(std::string("Unable to read alignment from \"")
                         + fn + std::string("\""));
  }
  if (i != tipCount){
    throw PLLException("Some taxa are missing from FASTA file");
  }
  return std::move(pm);
}

void ParsedMatrix::clear() {
  for (auto h : headers) {
    std::free(static_cast<void *>(const_cast<char *>(h))); // allocated in pll_fasta_getnext
  }
  for (auto s : seqData) {
    std::free(static_cast<void *>(const_cast<char *>(s))); // allocated in pll_fasta_getnext
  }
}

} // namespace
