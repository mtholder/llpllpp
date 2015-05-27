#if ! defined(__LLPLLPLUSPLUS_MODEL_STORAGE_DESCRIPTION_HPP__)
#define __LLPLLPLUSPLUS_MODEL_STORAGE_DESCRIPTION_HPP__
// An argument to create partitions that bundles info about
//  how to allocate storage for the calculations.

#include "llpllpp/base_includes.hpp"

namespace pllpp {

enum class DataCharEncodings {
  BINARY_DATA_ENCODING,
  NUCLEOTIDE_DATA_ENCODING,
  AMINO_ACID_DATA_ENCODING
};

const unsigned * getMapForEncoding(DataCharEncodings d);

class ModelStorageDescription {
  bool parameterizedByFullRateMatrix = true;
  public:
  DataCharEncodings dataEncoding;
  unsigned int numStates;
  unsigned int numRateCats;
  ArchAttribEnum archAttributes;
  ModelStorageDescription(DataCharEncodings encoding,
                          unsigned int nStates,
                          unsigned int nRateCats,
                          ArchAttribEnum arch)
    :dataEncoding(encoding),
    numStates(nStates),
    numRateCats(nRateCats),
    archAttributes(arch) {
    assert((nStates == 2 && encoding == DataCharEncodings::BINARY_DATA_ENCODING)
           || (nStates == 4 && encoding == DataCharEncodings::NUCLEOTIDE_DATA_ENCODING)
           || (nStates == 20 && encoding == DataCharEncodings::AMINO_ACID_DATA_ENCODING));
  }
  const unsigned * getDataEncodingMap() const {
    return getMapForEncoding(dataEncoding);
  }
  bool isParameterizedByFullRateMatrix() const {
    return parameterizedByFullRateMatrix;
  }
};

class FixedModelRef {
  private:
    const char * name;
    const double * stateFreq;
    const double * exchangeabilies;
  public:
    FixedModelRef(const char * modelName, const double *stFreq, const double * exch) 
      :name(modelName),
      stateFreq(stFreq),
      exchangeabilies(exch) {
    }
    const char * getName() const {
      return name;
    }
    const double * getStateFrequenciesPtr() const {
      return stateFreq;
    }
    const double * getExchangeabilityParamsPtr() const {
      return exchangeabilies;
    }
};
} // namespace pllpp
#endif
