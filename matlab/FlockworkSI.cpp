#include "CastResult.h"
#include "SI.h"
#include "ResultClasses.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t N, number_of_vaccinated, number_of_infected, seed;
    double Q, inf_rate, rew_rate;
    bool use_random_rew, eq_flw;
    
    if (nrhs!=10)
    {
        mexPrintf("Got %d input arguments. ", nrhs);
        throw length_error("Invalid number of input arguments. Has to be 10.");
    }

    if (nlhs!=4)
    {
        mexPrintf("Got %d output arguments. ", nrhs);
        throw length_error("Invalid number of output arguments. Has to be 4.");
    }


    vector < pair <size_t,size_t > > edgelist = get_edgelist(prhs[0]);

    read_single_value(prhs[1],N);
    read_single_value(prhs[2],Q);
    read_single_value(prhs[3],inf_rate);
    read_single_value(prhs[4],rew_rate);
    read_single_value(prhs[5],number_of_vaccinated);
    read_single_value(prhs[6],number_of_infected);
    read_single_value(prhs[7],use_random_rew);
    read_single_value(prhs[8],eq_flw);
    read_single_value(prhs[9],seed);

    /*
    if (eq_flw)
        mexPrintf("    eq_flw  =  true\n");
    else
        mexPrintf("    eq_flw  =  false\n");
    */

    SIS_result result = SI( edgelist,
                            N,
                            Q,
                            inf_rate,
                            rew_rate,
                            number_of_vaccinated,
                            number_of_infected,
                            use_random_rew,
                            eq_flw,
                            seed
                           );

    /*
    for(auto edge: edgelist)
        mexPrintf("%d %d\n", edge.first, edge.second);

    vector < pair <size_t,double > > blub;
    blub.push_back(make_pair(1,2));
    blub.push_back(make_pair(3,4));
    blub.push_back(make_pair(5,6));

    mxArray * result;
    result = cast_EPI_array_to_matlab(blub.begin(), blub.end());
    plhs[0] = result;
    */
    plhs[0] = cast_EPI_array_to_matlab(result.I_of_t.begin(), result.I_of_t.end());
    plhs[1] = cast_EPI_array_to_matlab(result.SI_of_t.begin(), result.SI_of_t.end());
    plhs[2] = cast_EPI_array_to_matlab(result.R0_of_t.begin(), result.R0_of_t.end());
    plhs[3] = cast_EPI_array_to_matlab(result.edge_list.begin(), result.edge_list.end());
}
