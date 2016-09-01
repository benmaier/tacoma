#include "SIS.h"
#include "ResultClasses.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"
#include "CastResult.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *rates_in_m, *result;
    const mwSize *dims;

    int dimx,dimy,numdims;

    rates_in_m = mxDuplicateArray(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimx = int( dims[0] );
    dimy = int( dims[1] );
    mexPrintf("Array has %d dimensions with N_1 = %d, N_2 = %d, entries.\n",numdims,dimx,dimy);

    vector <double> rates(dimx*dimy);


    vector < pair <size_t,double > > blub;
}
