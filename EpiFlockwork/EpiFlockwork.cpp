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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SIR.h"
#include "SIRS.h"
#include "SIS.h"
#include "ResultClasses.h"

using namespace std;
namespace py = pybind11;

PYBIND11_PLUGIN(EpiFlockwork) {
    py::module m("EpiFlockwork", "Module to equilibrate a flockwork in a fast manner");
    
    m.def("SIS", &SIS, "Simulate an SIS process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
            );

    m.def("SIR", &SIR, "Simulate an SIR process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("t_run_total") = 0,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
            );

    m.def("SIRS", &SIRS, "Simulate an SIRS process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("susceptible_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("t_run_total") = 0,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
            );

    py::class_<SIR_result>(m,"SIR_result")
        .def(py::init<>())
        .def_readwrite("I_of_t", &SIR_result::I_of_t)
        .def_readwrite("R_of_t", &SIR_result::R_of_t)
        .def_readwrite("SI_of_t", &SIR_result::SI_of_t)
        .def_readwrite("R0_of_t", &SIR_result::R0_of_t)
        .def_readwrite("edge_list", &SIR_result::edge_list);

    py::class_<SIS_result>(m,"SIS_result")
        .def(py::init<>())
        .def_readwrite("I_of_t", &SIS_result::I_of_t)
        .def_readwrite("SI_of_t", &SIS_result::SI_of_t)
        .def_readwrite("R0_of_t", &SIS_result::R0_of_t)
        .def_readwrite("edge_list", &SIS_result::edge_list);


    return m.ptr();

}
