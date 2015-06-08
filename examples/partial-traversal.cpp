#include "llpllpp/llpllpp.hpp"
#include <iostream>
#include <iomanip>
#include <cstdarg>
#include <ctime>
#include <pll.h>
using namespace pllpp;
int doMain(int argc, char * argv[]);
void calcLikeDemo(const std::string & newickFilename, const std::string &fastaFilename);

constexpr unsigned NUM_STATES = 4U;
constexpr unsigned NUM_RATE_CATS = 4U;

void calcLikeDemo(const std::string & newickFilename, const std::string &fastaFilename) {
  std::shared_ptr<UTree> tree{std::move(UTree::parseNewick(newickFilename))};
  tree->setMissingBranchLength(0.000001);
  auto inpMatrix = ParsedMatrix::parseFasta(fastaFilename, tree->getOTUSet());
  ModelStorageDescription msd{DataCharEncodings::NUCLEOTIDE_DATA_ENCODING,
                              NUM_STATES,
                              NUM_RATE_CATS,
                              ArchAttribEnum::LLPLL_ATTRIB_ARCH_SSE};
  UPhyloCalculator phyCalc(*inpMatrix, msd, tree);
  inpMatrix.release(); // we have copied the data into the phyCalc, and are done with the raw copy.
  // we only have one block of data, with index= 0
  auto & model = phyCalc.getModel(0);
  // initialize the array of base frequencies
  const std::vector<double> stateFreqs{ 0.17, 0.19, 0.25, 0.39};
  // substitution rates for the 4x4 GTR model. This means we need exactly
  //   (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal 
  const std::vector<double> subsParams{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  model.setStateFrequencies(stateFreqs);
  model.setExchangeabilityParams(subsParams);
  model.getRateHet().setAlphaOfGammaDist(1.0);
  constexpr unsigned int N_DRAWS = 10;
  const auto innerNodeArray= phyCalc.getInnerNodesArray();
  for (auto i = 0U; i < N_DRAWS; ++i) {
    int r = rand() % phyCalc.getInnerNodesCount();
    auto randNode = innerNodeArray[r];
    int randDir = rand() % 3;
    for (auto j = 0; j < randDir; ++j) {
      randNode = randNode->next;
    }
    phyCalc.setVirtualRoot(randNode);
    phyCalc.updateProbMatrices(0);
    phyCalc.updatePartials(0);
    std::cout << "Log-L: " << std::setprecision(9) << phyCalc.computeLogLikelihood(0) << '\n';
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
  std::srand(std::time(nullptr));
  return main_wrapper(argc, argv, doMain);
}
