mex -setup C++;
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./EpiFlockwork' ./matlab/FlockworkSIS.cpp;