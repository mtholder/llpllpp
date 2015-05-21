#if ! defined(__LLPLLPLUSPLUS_PLL_EXCEPTION_HPP__)
#define __LLPLLPLUSPLUS_PLL_EXCEPTION_HPP__
#include <exception>
#include <string>
#include "llpllpp/base_includes.hpp"
namespace pllpp {

class PLLException: public std::exception {
  const std::string msg;
  public:
  PLLException(const std::string m)
    :msg(m) {
  }
  virtual const char* what() const noexcept {
    return msg.c_str();
  }
};

} // namespace 
#endif
