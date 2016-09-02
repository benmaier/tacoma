mex -setup C++;
%mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./EpiFlockwork -I./matlab' matlab/FlockworkSIS.cpp matlab/CastResult.cpp EpiFlockwork/SIS.cpp EpiFlockwork/Utilities.cpp EpiFlockwork/Events.cpp;
%mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./EpiFlockwork -I./matlab' matlab/FlockworkSIRS.cpp matlab/CastResult.cpp EpiFlockwork/SIRS.cpp EpiFlockwork/Utilities.cpp EpiFlockwork/Events.cpp;
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./EpiFlockwork -I./matlab' matlab/FlockworkSIR.cpp matlab/CastResult.cpp EpiFlockwork/SIR.cpp EpiFlockwork/Utilities.cpp EpiFlockwork/Events.cpp;
