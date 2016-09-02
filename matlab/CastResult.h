/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */
#ifndef __CASTRESULT_H__
#define __CASTRESULT_H__

#include "math.h"
#include "matrix.h"
#include "mex.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <random>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <cstdlib>
#include <tuple>


using namespace std;

vector < pair<size_t,size_t> > get_edgelist(const mxArray *m_edges);

template < class fwditer >
mxArray * cast_EPI_array_to_matlab( fwditer begin, fwditer end )
{
    fwditer current_pair = begin;
    size_t len_container = distance(begin, end);
    mxArray *result;
    result = mxCreateDoubleMatrix(len_container,2,mxREAL);
    double *r;
    r = mxGetPr(result);

    size_t index = 0;
    while (current_pair != end) 
    {
        r[index] = (double) (*current_pair).first;
        r[index+len_container] = (double) (*current_pair).second;

        current_pair++;
        index++;
    }

    return result;
}

template < typename num >
void read_single_value(const mxArray *m_single, num &val)
{
    double *value;
    size_t numdims, dimx, dimy;
    const mwSize *dims;

    value = mxGetPr(m_single);
    dims = mxGetDimensions(m_single);
    numdims = mxGetNumberOfDimensions(m_single);
    dimx = int( dims[0] );
    dimy = int( dims[1] );
    
    if ( (dimx!=1) || (dimy!=1) )
        throw length_error("Got a non-single value argument as a single value argument.");

    val = (num) *value;

}

#endif
