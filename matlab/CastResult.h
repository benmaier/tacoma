#include "math.h"
#include "matrix.h"
#include "mex.h"
using namespace std;

template < class fwditer >
mxArray * cast_EPI_array_to_matlab( fwditer begin, fwditer end );

vector < pair<size_t,size_t> > get_edgelist(mxArray &m_edges);
