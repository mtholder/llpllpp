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
set -x
CPPFLAGS="-I${PLL_ROOT_INC}" \
  CXX=$(which clang++) \
  CC=$(which clang) \
  LDFLAGS="-L${PLL_ROOT_LIB} -lpll -lm" \
  CXXFLAGS="-Wno-c++98-compat -Weverything -Wpadded -pedantic -g -O0 -std=c++1y -cxx-isystem=/usr/lib/gcc/x86_64-linux-gnu/4.9" \
  ../configure --prefix=$PWD/installed
