#include "llpllpp/aa_models.hpp"
#include "pll.h"
#include <cstdlib>
namespace pllpp {

std::vector<FixedModelRef> AAModel::fixedModels;
/* from libpll/examples/protein-list/pl */
constexpr unsigned PROT_MODELS_COUNT = 19;

static const double * protein_models_rates_list[PROT_MODELS_COUNT] = {
   pll_aa_rates_dayhoff,
   pll_aa_rates_lg,
   pll_aa_rates_dcmut,
   pll_aa_rates_jtt,
   pll_aa_rates_mtrev,
   pll_aa_rates_wag,
   pll_aa_rates_rtrev,
   pll_aa_rates_cprev,
   pll_aa_rates_vt,
   pll_aa_rates_blosum62,
   pll_aa_rates_mtmam,
   pll_aa_rates_mtart,
   pll_aa_rates_mtzoa,
   pll_aa_rates_pmb,
   pll_aa_rates_hivb,
   pll_aa_rates_hivw,
   pll_aa_rates_jttdcmut,
   pll_aa_rates_flu,
   pll_aa_rates_stmtrev
 };

static const double * protein_models_freqs_list[PROT_MODELS_COUNT] = {
   pll_aa_freqs_dayhoff,
   pll_aa_freqs_lg,
   pll_aa_freqs_dcmut,
   pll_aa_freqs_jtt,
   pll_aa_freqs_mtrev,
   pll_aa_freqs_wag,
   pll_aa_freqs_rtrev,
   pll_aa_freqs_cprev,
   pll_aa_freqs_vt,
   pll_aa_freqs_blosum62,
   pll_aa_freqs_mtmam,
   pll_aa_freqs_mtart,
   pll_aa_freqs_mtzoa,
   pll_aa_freqs_pmb,
   pll_aa_freqs_hivb,
   pll_aa_freqs_hivw,
   pll_aa_freqs_jttdcmut,
   pll_aa_freqs_flu,
   pll_aa_freqs_stmtrev
 };

static const char * protein_models_names_list[PROT_MODELS_COUNT] = {
   "DAYHOFF",
   "LG",
   "DCMUT",
   "JTT",
   "MTREV",
   "WAG",
   "RTREV",
   "CPREV",
   "VT",
   "BLOSUM62",
   "MTMAM",
   "MTART",
   "MTZOA",
   "PMB",
   "HIVB",
   "HIVW",
   "JTTDCMUT",
   "FLU",
   "STMTREV"
};

const std::vector<FixedModelRef> & AAModel::getFixedModels() {
  if (fixedModels.empty()) {
    fixedModels.reserve(PROT_MODELS_COUNT);
    for (auto i = 0U; i < PROT_MODELS_COUNT; ++i) {
      fixedModels.emplace_back(protein_models_names_list[i],
                               protein_models_freqs_list[i],
                               protein_models_rates_list[i]);
    }
  }
  return fixedModels;
}

} // namespace

