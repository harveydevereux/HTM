#!/bin/bash
g++ -std=c++11 tests.cpp -o tests && ./tests --success
g++ -O3 -std=c++14 HTM.cpp -o HTM && ./HTM
julia-1.3 plots.jl
