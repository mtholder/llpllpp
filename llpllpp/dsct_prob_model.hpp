#if ! defined(__LLPLLPLUSPLUS_DSCTPROB_MODEL_HPP__)
#define __LLPLLPLUSPLUS_DSCTPROB_MODEL_HPP__
// A discrete-state, continuous time probability model calculator
#include <vector>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/model_storage_description.hpp"
namespace pllpp {

class RateHetModel {
  mutable std::vector<double> rates;
  double alphaParam;
  unsigned long rateSetCounter = 0;
  mutable unsigned long rateCalcCounter = 0;
  public:
  RateHetModel(std::size_t numCategories)
    :rates(numCategories, -1.0) {
  }
  void setAlphaOfGammaDist(double alpha) {
    rateSetCounter += 1;
    assert(alpha > 0.0);
    alphaParam = alpha;
  }
  unsigned long getCounter() const noexcept {
    return rateSetCounter;
  }
  const std::vector<double> & getRates() const;
};

class DSCTProbModel {
  RateHetModel rateHet;
  std::vector<double> stateFrequencies;
  std::vector<double> exchangeParameters;
  unsigned long stateFreqSetCounter = 0;
  unsigned long exchangeSetCounter = 0;
  
  public:
  DSCTProbModel(const ModelStorageDescription & msd)
    :rateHet(msd.numRateCats),
    stateFrequencies(msd.numStates, 1.0/(static_cast<double>(msd.numStates))) {
    stateFreqSetCounter = 1;
  }
  DSCTProbModel(const ModelStorageDescription & msd,
                const std::vector<double> &stateFreqVec,
                const std::vector<double> &exchangeParams)
    :rateHet(msd.numRateCats),
    stateFrequencies(stateFreqVec),
    exchangeParameters(exchangeParams) {
    stateFreqSetCounter = 1;
    exchangeSetCounter = 1;
    if (stateFrequencies.size() != msd.numStates) {
      throw PLLException("Number of element in stateFreqVec != numStates");
    }
  }
  unsigned long getStateFreqCounter() const noexcept {
    return stateFreqSetCounter;
  }
  unsigned long getExchangeCounter() const noexcept {
    return exchangeSetCounter;
  }
  const RateHetModel & getRateHet() const {
    return rateHet;
  }
  RateHetModel & getRateHet() {
    return rateHet;
  }
  void setStateFrequencies(const std::vector<double> v) {
    assert(v.size() == stateFrequencies.size());
    stateFreqSetCounter += 1;
    stateFrequencies = v;
  }
  void setExchangeabilityParams(const std::vector<double> v) {
    assert(exchangeParameters.empty() || v.size() == exchangeParameters.size());
    exchangeSetCounter += 1;
    exchangeParameters = v;
  }
  const std::vector<double> & getStateFrequencies() const {
    return stateFrequencies;
  }
  const std::vector<double> & getExchangeabilityParams() const {
    return exchangeParameters;
  }
};

} // namespace pllpp
#endif
