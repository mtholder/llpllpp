#include "llpllpp/phylo_calculator.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {
const std::vector<double> & RateHetModel::getRates() const {
  if (rateCalcCounter != rateSetCounter) {
    pll_compute_gamma_cats(alphaParam, static_cast<int>(rates.size()), &rates[0]);
    rateCalcCounter = rateSetCounter;
  }
  return rates;
}

} // namespace

