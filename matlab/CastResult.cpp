#include "CastResult.h"
#include <stdexcept>

using namespace std;

vector < pair<size_t,size_t> > get_edgelist(const mxArray * m_edges)
{
    double *m_edgelist;
    size_t numdims, dimx, dimy;
    const mwSize *dims;

    m_edgelist = mxGetPr(m_edges);
    dims = mxGetDimensions(m_edges);
    numdims = mxGetNumberOfDimensions(m_edges);
    dimx = int( dims[0] );
    dimy = int( dims[1] );

    //if (numdims != 2)
    //{
    //    mexPrintf("Wrong dimensions in edgelist: [%d,%d]",dimx,dimy);
    //    throw length_error();
    //}

    vector < pair<size_t,size_t> > edges;

    if (dimx>0)
    {
        if (dimy != 2)
        {
            mexPrintf("Wrong dimensions in edgelist: [%d,%d]\n",dimx,dimy);
            throw length_error("");
        }

        for(size_t x=0; x<dimx; x++)
        {
            size_t p[dimy];
 
            for(size_t y=0; y<dimy; y++)
                p[y] = m_edgelist[x+y*dimx];

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
        mexPrintf("Wrong dimensions in edgelist: [%d,%d]\n",dimx,dimy);
        throw length_error("");
    }

}
