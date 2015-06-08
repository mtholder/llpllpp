#include "llpllpp/llpllpp.hpp"
#include <iostream>
#include <iomanip>
using namespace pllpp;
int doMain(int argc, char * argv[]);
void calcLikeDemo(const std::string & newickFilename, const std::string &fastaFilename);

constexpr unsigned NUM_STATES = 20U;
constexpr unsigned NUM_RATE_CATS = 4U;


void calcLikeDemo(const std::string & newickFilename, const std::string &fastaFilename) {
  std::shared_ptr<UTree> tree{std::move(UTree::parseNewick(newickFilename))};
  tree->setMissingBranchLength(0.000001);
  auto inpMatrix = ParsedMatrix::parseFasta(fastaFilename, tree->getOTUSet());
  ModelStorageDescription msd{DataCharEncodings::AMINO_ACID_DATA_ENCODING,
                              NUM_STATES,
                              NUM_RATE_CATS,
                              ArchAttribEnum::LLPLL_ATTRIB_ARCH_SSE};
  UPhyloCalculator phyCalc(*inpMatrix, msd, tree);
  inpMatrix.release(); // we have copied the data into the phyCalc, and are done with the raw copy.
  // we only have one block of data, with index= 0
  auto & model = phyCalc.getModel(0);
  model.getRateHet().setAlphaOfGammaDist(1.0);
  for (auto aaMod : AAModel::getFixedModels()) {
    model.setStateFrequenciesPtr(aaMod.getStateFrequenciesPtr());
    model.setExchangeabilityParamsPtr(aaMod.getExchangeabilityParamsPtr());
    phyCalc.updateProbMatrices(0);
    phyCalc.updatePartials(0);
    std::cout << "Log-L (" << aaMod.getName() << "): "
              << std::setprecision(9) << phyCalc.computeLogLikelihood(0) << '\n';
  }
}


int doMain(int argc, char * argv[]) {
  if (argc != 3) {
    throw PLLException("Expects 2 arguments: [newick] [fasta]");
  }
  const std::string newickFilename{argv[1]};
  const std::string fastaFilename{argv[2]};
  calcLikeDemo(newickFilename, fastaFilename);
  return 0;
}

int main(int argc, char * argv[]) {
  return main_wrapper(argc, argv, doMain);
}
