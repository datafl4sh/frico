/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025,2026 Matteo Cicuttin - IV3IWE
 * Politecnico di Torino
 * Dipartimento di Scienze Matematiche "G. L. Lagrange"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <complex>
#include <limits>

#include "emw_solver.h"
#include "sources.h"
namespace frico::maxwell {

struct port_values {
    std::complex<double>    voltage = 0.0;
    std::complex<double>    current = 0.0;
    std::complex<double>    power = 0.0;
    std::complex<double>    impedance = 0.0;
    std::complex<double>    gamma = 0.0;
    double                  swr = std::numeric_limits<double>::infinity();
};
struct antenna_data {
    std::vector<double>     gain_XY;
    std::vector<double>     gain_YZ;
    std::vector<double>     gain_XZ;
    port_values             port;
};
struct antenna_analysis {
    delta_gap   dgap;           /* delta-gap excitation */
    double      rdiag_dist;     /* distance for radiation diag computation */
    std::vector<antenna_data>   antdata;
};

bool init_from_spec(const char *, antenna_analysis& anta);
bool do_sweep(simulation& sim, antenna_analysis& ana);

} //namespace frico::maxwell