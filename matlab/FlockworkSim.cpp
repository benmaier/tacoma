#include "CastResult.h"
#include "SIS.h"
#include "SIR.h"
#include "SIRS.h"
#include "EqFlockwork.h"
#include "ResultClasses.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t N, t_run_total, seed;
    double Q;
    
    if (nrhs!=5)
    {
        mexPrintf("Got %d input arguments. ", nrhs);
        throw length_error("Invalid number of input arguments. Has to be 5.");
    }

    if (nlhs!=1)
    {
        mexPrintf("Got %d output arguments. ", nrhs);
        throw length_error("Invalid number of output arguments. Single output variable expected.");
    }

    vector < pair <size_t,size_t> > edgelist = get_edgelist(prhs[0]);

    read_single_value(prhs[1],N);
    read_single_value(prhs[2],Q);
    read_single_value(prhs[3],seed);
    read_single_value(prhs[4],t_run_total);

    vector < pair <size_t,size_t> > new_edgelist;
    new_edgelist = simulate_flockwork(edgelist,
                                      N,
                                      Q,
                                      seed,
                                      t_run_total
                                     );

    plhs[0] = cast_EPI_array_to_matlab(new_edgelist.begin(), new_edgelist.end());
}
