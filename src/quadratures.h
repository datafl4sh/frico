/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025, Matteo Cicuttin - IV3IWE
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

#include <vector>
#include "geom_point.h"
#include "geom_mesh.h"

namespace frico {

struct quadrature_point {
    point   p;
    double  w;
};

std::vector<quadrature_point>
dunavant(size_t, const point&, const point&, const point&);

std::vector<quadrature_point>
integrate(const mesh&, const triangle&, size_t degree);

std::vector<quadrature_point>
integrate_subtri(const mesh&, const triangle&, size_t degree);

} // namespace frico
