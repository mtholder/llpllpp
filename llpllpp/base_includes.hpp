#if !defined __LLPLLPLUSPLUS_BASE_INCLUDES_HPP__
#define __LLPLLPLUSPLUS_BASE_INCLUDES_HPP__
// This file is included by any header in the library.
//  It holds the very few defines that we use and  includes the logger. 

#include <cassert>
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored  "-Wweak-vtables"
#define NOT_IMPLEMENTED assert("not implemented"[0] == 'f');
#define UNREACHABLE assert(false);
#include "llpllpp/easylogging++.hpp"

//
enum class ArchAttribEnum {
  LLPLL_ATTRIB_ARCH_SSE =    0x01,
  LLPLL_ATTRIB_ARCH_AVX =    0x02,
  LLPLL_ATTRIB_ARCH_AVX2 =   0x04,
  LLPLL_ATTRIB_ARCH_AVX512 = 0x08
};

#include "llpllpp/pll_exception.hpp"


#endif
