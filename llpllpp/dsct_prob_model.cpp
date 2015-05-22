#include "llpllpp/phylo_calculator.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {
const std::vector<double> & RateHetModel::getRates() const {
  if (rateCalcCounter != rateSetCounter) {
    double * r = &(rates[0]);
    pll_compute_gamma_cats(alphaParam, static_cast<int>(rates.size()), r);
    rateCalcCounter = rateSetCounter;
  }
  return rates;
}

} // namespace

