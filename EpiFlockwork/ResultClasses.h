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

#ifndef __RESULTCLASSES_H__
#define __RESULTCLASSES_H__

#include "Events.h"
#include "Utilities.h"
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
#include <tuple>

struct SIR_result 
{
    SIR_result(
                const vector < pair < double, size_t > > & I_of_t,  
                const vector < pair < double, size_t > > & R_of_t, 
                const vector < pair < double, size_t > > & SI_of_t, 
                const vector < pair < double, double > > & R0_of_t,
                const set < pair < size_t, size_t >  > & edge_list
            )
    {
       I_of_t(I_of_t); 
       R_of_t(R_of_t); 
       SI_of_t(SI_of_t); 
       R0_of_t(R0_of_t); 
       edge_list(edge_list);
    }

    const vector < pair < double, size_t > > &I_of_t const { return I_of_t; }
    const vector < pair < double, size_t > > &R_of_t const { return R_of_t; }
    const vector < pair < double, size_t > > &SI_of_t const { return SI_of_t; }
    const vector < pair < double, double > > &R0_of_t const { return R0_of_t; }
    const set < pair < size_t, size_t >  > &edge_list() const { return edge_list; }

    vector < pair < double, size_t > > I_of_t; 
    vector < pair < double, size_t > > R_of_t; 
    vector < pair < double, size_t > > SI_of_t; 
    vector < pair < double, double > > R0_of_t;
    set < pair < size_t, size_t >  > edge_list;
    
}

#endif
