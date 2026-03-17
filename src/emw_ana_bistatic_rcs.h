/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025-2026, Matteo Cicuttin - IV3IWE
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

#include <fstream>
#include "emw_solver.h"
namespace frico::maxwell {
struct bistatic_rcs_analysis {
    double          radar_R = 10.0;
    double          radar_theta = 0.0;
    double          radar_phi = 0.0;
    double          Etheta = 1.0;
    double          Ephi = 1.0;
    std::ofstream   ofs;
};

bool init_from_spec(const char *spec, bistatic_rcs_analysis& ana);

bool do_sweep(simulation& sim, const bistatic_rcs_analysis& ana);

} // namespace frico::maxwell