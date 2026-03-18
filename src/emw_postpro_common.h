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

#include <expected>

#include "eigen.h"
#include "emw_solver.h"
#include "utils.h"

namespace frico::maxwell {

std::pair<ezvec3, ezvec3> eval_fields(const simulation&, size_t, const point&);
void eval_fields(const simulation&, size_t, const mesh&, zdfield&, zdfield&);
std::complex<double> gamma(std::complex<double>, double);
double swr(std::complex<double>, double);
double swr(std::complex<double>);
bool make_sampling_grid(frico::mesh& msh, const frico::point& c,
    double r, double h);
void write_fields(simulation& sim, size_t ctx_number);


} // namespace frico::maxwell

