#if ! defined(__LLPLLPLUSPLUS_DSCTPROB_MODEL_HPP__)
#define __LLPLLPLUSPLUS_DSCTPROB_MODEL_HPP__
// A discrete-state, continuous time probability model calculator
#include <vector>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/model_storage_description.hpp"
namespace pllpp {

class RateHetModel {
  std::vector<double> rates;
  double alphaParam;
  public:
  RateHetModel(std::size_t numCategories)
    :rates(numCategories, -1.0) {
  }
  void setAlphaOfGammaDist(double alpha) {
    assert(alpha > 0.0);
    alphaParam = alpha;
  }

};
class DSCTProbModel {
  RateHetModel rateHet;
  std::vector<double> stateFrequencies;
  std::vector<double> exchangeParameters;

  public:
  DSCTProbModel(const ModelStorageDescription & msd)
    :rateHet(msd.numRateCats),
    stateFrequencies(msd.numStates, 1.0/(static_cast<double>(msd.numStates))) {
  }
  DSCTProbModel(const ModelStorageDescription & msd,
                const std::vector<double> &stateFreqVec,
                const std::vector<double> &exchangeParams)
    :rateHet(msd.numRateCats),
    stateFrequencies(stateFreqVec),
    exchangeParameters(exchangeParams) {
    if (stateFrequencies.size() != msd.numStates) {
      throw PLLException("Number of element in stateFreqVec != numStates");
    }
  }
  RateHetModel & getRateHet() {
    return rateHet;
  }
  void setStateFrequencies(const std::vector<double> v) {
    assert(v.size() == stateFrequencies.size());
    stateFrequencies = v;
  }
  void setExchangeabilityParams(const std::vector<double> v) {
    assert(exchangeParameters.empty() || v.size() == exchangeParameters.size());
    exchangeParameters = v;
  }
};

} // namespace pllpp
#endif
