# llpllpp

Note: libpll for the code base that provides  the real functionality!
That lib is at: http://www.libpll.org/

However, this repo is wrapping an as-yet-unreleased low-level subset of PLL.
This is just a C++ wrapper around the low-level PLL.

Currently this relies on the "minor" branch of that repo (but only because of
some `const` stuff - it is not far from working on the master branch).

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

You should get:

    ./examples/newick-fasta-unrooted ../data/small.tre ../data/small.fas
    Log-L: -5895.90491

