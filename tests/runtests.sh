#!/bin/bash
BOLD=$(tput bold)
normal=$(tput sgr0)
GREEN='\033[0;32m'
NORMAL='\033[0m'

echo -e "${GREEN}${BOLD}Compiling and running Catch.hpp test suite...${normal}${NORMAL}"
g++ -std=c++11 tests.cpp -o tests && ./tests --success
echo -e "${GREEN}${BOLD}Compiling and running benchmarks and examples...${normal}${NORMAL}"
g++ -O3 -std=c++14 HTM.cpp -o HTM && ./HTM
echo -e "${GREEN}${BOLD}Plotting benchmarks and examples in Julia...${normal}${NORMAL}"
julia-1.3 plots.jl
