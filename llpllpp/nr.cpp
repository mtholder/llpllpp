#include "llpllpp/nr.hpp"

#if 0
/* the function below actually implements the iterative Newton-Raphson procedure.
   It is particularly messy and hard to read because for the case of per-partition branch length 
   estimates it needs to keep track of whetehr the Newton Raphson procedure has 
   converged for each partition individually. 
   The rationale for doing it like this is also provided in:
   A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009,

*/
static void topLevelMakenewz(MTInstance & instance,
                             const std::vector<double> & z0, // init branch lengths
                             int _maxiter,
                             std::vector<double> & result) {
    ConcatModelInfo & si = instance.GetConcatModelInfo();
    const auto numBranchParts = z0.size();
      // initialize loop convergence variables etc. 
      // maxiter is the maximum number of NR iterations we are going to do before giving up 
    std::vector<double> z = z0;
    std::vector<int> maxiter(numBranchParts, _maxiter);
    std::deque<bool> outerConverged(numBranchParts, false);
    std::deque<bool> curvatOK(numBranchParts, true);
    std::vector<double> coreLogZ(numBranchParts);
    std::vector<double> zstep(numBranchParts);
    std::vector<double> zprev(numBranchParts);
    for (;;) { // check if we are done for partition i or if we need to adapt the branch length again */
        for(auto i = 0U; i < numBranchParts; i++) {
            if (!outerConverged[i]) {
                if (curvatOK[i]) {
                    curvatOK[i] = false;
                    zprev[i] = z[i];
                    zstep[i] = (1.0 - PLL_ZMAX) * z[i] + PLL_ZMIN;
                }
                    // other case, the outer loop hasn't converged but we are trying to approach the maximum from the wrong side 
                enforce_bound(z[i], PLL_ZMIN, PLL_ZMAX);
                coreLogZ[i] = std::log(z[i]);
            }
        SetExecuteMaskByBranchLengthPart(instance, curvatOK, true);
        StoreExecuteMaskInTraversalDescriptor(instance);
        StoreValuesInTraversalDescriptor(instance, coreLogZ);
        std::vector<double> dlnLdlz;
        std::vector<double> d2lnLdlz2;
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
            // if this is the first iteration of NR we will need to first do this one-time call 
            // of maknewzIterative() Note that, only this call requires broadcasting the traversal descriptor,
            // subsequent calls to pllMasterBarrier(PLL_THREAD_MAKENEWZ, tr); will not require this
        if(firstIteration) {
            tr->td[0].traversalHasChanged = true; 
            pllMasterBarrier (tr, pr, PLL_THREAD_MAKENEWZ_FIRST);
            firstIteration = false; 
            tr->td[0].traversalHasChanged = false; 
        } else {
            pllMasterBarrier(tr, pr, PLL_THREAD_MAKENEWZ);
        }
        branchLength_parallelReduce(tr, (double*)dlnLdlz, (double*)d2lnLdlz2, numBranchParts);
#else
            // sequential part, if this is the first newton-raphson implementation,
            // do the precomputations as well, otherwise just execute the computation
            // of the derivatives
        if (firstIteration) {
            makenewzIterative(tr, pr);
            firstIteration = false;
        }
        execCore(tr, pr, dlnLdlz, d2lnLdlz2);
#endif
            // do a NR step, if we are on the correct side of the maximum that's okay, otherwise 
            // shorten branch 
        for(i = 0; i < numBranchParts; i++) {
            if(!outerConverged[i] && !curvatOK[i]) {
                if ((d2lnLdlz2[i] >= 0.0) && (z[i] < PLL_ZMAX)) {
                    zprev[i] = z[i] = 0.37 * z[i] + 0.63; // Bad curvature, shorten branch
                } else {
                    curvatOK[i] = true;
                }
            }
        }
            // do the standard NR step to obrain the next value, depending on the state for eahc partition
        for(i = 0; i < numBranchParts; i++) {
            if (curvatOK[i] && !outerConverged[i]) {
                const double pb = 0.25 * zprev[i] + 0.75;
                if (d2lnLdlz2[i] < 0.0) {
                    double tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
                    if (tantmp < 100) {
                        z[i] *= exp(tantmp);
                    }
                    enforce_min_bound(z[i], PLL_ZMIN);
                    enforce_max_bound(z[i], pb);
                } else {
                    z[i] = pb;
                }
            }
            enforce_max_bound(z[i], PLL_ZMAX);
                // decrement the maximum number of itarations */
            maxiter[i] = maxiter[i] - 1;
                // check if the outer loop has converged */
            if(std::abs(z[i] - zprev[i]) > zstep[i]) {
                   // We should make a more informed decision here, based on the log like improvement
                if (maxiter[i] < -20) {
                    z[i] = z0[i];
                    outerConverged[i] = true;
                } else {
                    outerConverged[i] = false;
                }
            } else {
                outerConverged[i] = true;
            }
        }
        if (all_true(outerConverged)) {
            break;
        }
    }
    SetExecuteMaskForAllPart(instance, true);
    result = z;
}
#endif
