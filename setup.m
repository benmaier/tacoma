% use C++ in the following
mex -setup C++;

% compile functions
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./cFlockwork -I./matlab' matlab/FlockworkSIS.cpp matlab/CastResult.cpp cFlockwork/SIS.cpp cFlockwork/Utilities.cpp cFlockwork/Events.cpp cFlockwork/EqFlockwork.cpp;
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./cFlockwork -I./matlab' matlab/FlockworkSIRS.cpp matlab/CastResult.cpp cFlockwork/SIRS.cpp cFlockwork/Utilities.cpp cFlockwork/Events.cpp cFlockwork/EqFlockwork.cpp;
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./cFlockwork -I./matlab' matlab/FlockworkSIR.cpp matlab/CastResult.cpp cFlockwork/SIR.cpp cFlockwork/Utilities.cpp cFlockwork/Events.cpp cFlockwork/EqFlockwork.cpp;
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./cFlockwork -I./matlab' matlab/FlockworkSim.cpp matlab/CastResult.cpp cFlockwork/SIR.cpp cFlockwork/Utilities.cpp cFlockwork/Events.cpp cFlockwork/EqFlockwork.cpp;
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./cFlockwork -I./matlab' matlab/FlockworkEq.cpp matlab/CastResult.cpp cFlockwork/SIR.cpp cFlockwork/Utilities.cpp cFlockwork/Events.cpp cFlockwork/EqFlockwork.cpp;

% move compiled functions to new folder
mkdir matlabbuild
movefile('FlockworkS*','./matlabbuild');
movefile('FlockworkEq*','./matlabbuild');

% add path to matlab environment via startup-file in user directory
up = userpath;
startuppath = [up(1:end-1),'/startup.m'];
libpath = [pwd,'/matlabbuild'];

fid = fopen(startuppath, 'at');  % append to possibly existing startup-file
fprintf(fid,['\naddpath(''',libpath,''')\n']);
fclose(fid);

% add path for this session
addpath(libpath);