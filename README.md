# llpllpp
Is a C++ wrapper around the low-level PLL lib

See https://github.com/xflouris/libpll for the code base that provides 
  the real functionality!
[![Build Status](https://magnum.travis-ci.com/xflouris/libpll.svg?token=rjft2y6GBHow4SDyjuoy&branch=master)](https://magnum.travis-ci.com/xflouris/libpll)
[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# Installation

## Prerequisites

  1. build the low-level pll from https://github.com/xflouris/libpll
  2. Set `PLL_ROOT_INC` and `PLL_ROOT_LIB` in your env to point to the `src` subdirectory of `libpll`.
  3. Add `PLL_ROOT_LIB` to your `LD_LIBRARY_PATH` (or equivalent)

You also need a modern C++ compiler and autotools

## One time-bootstrapping

  1. run `sh bootstrap.sh` to create the `configure` script
  2. `mkdir build; cd build` if you want to keep build products out of your source directories 
  3. run configure. You may be able to get away with `bash ../reconf-clang.sh` or `bash ../reconf-gcc.sh`
      consult those scripts to see how various env variables and arguments are passed to `configure`
      and tweak for your system.

## building

    make
    make check
    make install
    make installcheck

(though there currently aren't any tests)