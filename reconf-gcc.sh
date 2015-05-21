#!/bin/bash
if test -z $PLL_ROOT_INC
then
    echo PLL_ROOT_INC must be in your env
    exit 1
fi
if test -z $PLL_ROOT_LIB
then
    echo PLL_ROOT_LIB must be in your env
    exit 1
fi
CPPFLAGS="-I${PLL_ROOT_INC}" \
  CXX=$(which g++) \
  CC=$(which gcc) \
  LDFLAGS="-L${PLL_ROOT_LIB} -lpll -lm" \
  CXXFLAGS="-Waddress -Warray-bounds -Wc++11-compat -Wchar-subscripts -Wcomment -Wformat -Wmain -Wmaybe-uninitialized -Wmissing-braces -Wnonnull -Wopenmp-simd -Wparentheses -Wreorder -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow -Wswitch -Wtrigraphs -Wuninitialized -Wno-unknown-pragmas -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var -Wno-pragmas -pedantic -g -O0 -std=c++14" \
    ../configure --prefix=$PWD/installed
