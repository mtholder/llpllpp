#if ! defined(__LLPLLPLUSPLUS_UTIL_HPP__)
#define __LLPLLPLUSPLUS_UTIL_HPP__
#include <functional>
#include <string>
#include <iostream>
#include "llpllpp/base_includes.hpp"
#include "llpllpp/pll_exception.hpp"

namespace pllpp {
// Wraps main with a catch to report the exception before returning 1
int main_wrapper(int argc, char * argv[], std::function<int(int, char * [])> fn);



// Wraps main with a catch to report the exception before returning 1
inline int main_wrapper(int argc, char * argv[], std::function<int(int, char * [])> fn) {
  std::string program;
  if (argc > 0) {
    program.assign(argv[0]);
  }
  try {
    return fn(argc, argv);
  } catch (const PLLException & x) {
    std::cerr << program << ": " << x.what() << std::endl;
    return 1;
  } catch (const std::exception & y) {
    std::cerr << program << ": Unexpected exception thrown. " << y.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << program << ": Unknown exception type thrown." << std::endl;
    return 1;
  }
}

}// namespace
#endif
