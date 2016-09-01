#include "CastResult.h"

using namespace std;

template < class fwditer >
mxArray * cast_EPI_array_to_matlab( fwditer begin, fwditer end )
{
    fwditer current_pair = begin;
    size_t len_container = distance(begin, end);
    mxArray *result;
    result = mxCreateDoubleMatrix(2,len_container,mxREAL);
    double *r;
    r = mxGetPr(result);

    while (current_pair != end) 
    {
        *r = (double) (*current_pair).first;
        r++;

        *r = (double) (*current_pair).second;
        r++;

        //advance(current_pair,1);
        current_pair++;
    }

    return result;
}

vector < pair<size_t,size_t> > get_edgelist(mxArray &m_edges)
{
    double *m_edgelist;
    m_edgelist = mxGetPr(m_edges);
    dims = mxGetDimensions(m_edges);
    numdims = mxGetNumberOfDimensions(m_edges);
    dimx = int( dims[0] );
    dimy = int( dims[1] );

    vector < tuple<size_t,size_t> > edges;

    if (dimx>0)
    {
        if (dimy != 2)
        {
            mexPrintf("Wrong dimensions in edgelist: [%d,%d]",dimx,dimy);
            throw length_error();
        }

        for(size_t x=0; x<dimx; x++)
        {
            size_t p[dimy];
 
            for(size_t y=0; y<dimy; y++)
                p[y] = m_edgelist[y+x*dimy]

            edges.push_back( make_pair(p[0], p[1]) );
        }

        return edges;
    }
    else if (dimy == 0)
    {
        return edges;
    }
    else
    {
        mexPrintf("Wrong dimensions in edgelist: [%d,%d]",dimx,dimy);
        throw length_error();
    }

}
