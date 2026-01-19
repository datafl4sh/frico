/*
 * MoMmy - My experimental Method of Moments code
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

#include <iostream>
#include <optional>
#include <vector>
#include <array>

#include "geom_mesh.h"
#include "geom_point.h"
#include "eigen.hpp"

namespace mommy {

/* The convention used for the basis functions is the one of 
 * "MOM3D Method of Moments Code Theory Manual" by J. F. Shaeffer.
 * Rho+ and Rho- always point from the node to the shared edge,
 * T+ has sign +1 and T- has sign -1.
 * Beware that Gibson uses a different, opposite convention.
 */
struct basis_function
{
    double  Aminus;         // area of T-
    double  Aplus;          // area of T+
    size_t  itminus;        // index of T-
    size_t  itplus;         // index of T+
    point   pminus;         // vertex of T- opposite to the edge
    point   pplus;          // vertex of T+ opposite to the edge
    double  length;         // length of the edge
    size_t  edge_index;     // global edge index
    size_t  matrix_index;   // index of the edge in the matrix
    std::optional<size_t>   interface; // if internal interface, its tag

    vec3 rho_minus(const point& r) const;
    vec3 rho_plus(const point& r) const;

    edvec3 eval_minus(const point& r) const;
    edvec3 eval_plus(const point& r) const;

    double div_minus(const point& r) const;
    double div_plus(const point& r) const;
};

inline vec3
basis_function::rho_minus(const point& r) const
{
    return r - pminus;
}

inline vec3
basis_function::rho_plus(const point& r) const
{
    return r - pplus;
}

inline edvec3
basis_function::eval_minus(const point& r) const
{
    auto f = -length*rho_minus(r)/(2.0*Aminus);
    return { f.x(), f.y(), f.z() };
}

inline edvec3
basis_function::eval_plus(const point& r) const
{
    auto f = length*rho_plus(r)/(2.0*Aplus);
    return { f.x(), f.y(), f.z() };
}

inline double
basis_function::div_minus(const point& r) const
{
    return -length/(2.0*Aminus);
}

inline double
basis_function::div_plus(const point& r) const
{
    return length/(2.0*Aplus);
}

inline std::ostream&
operator<<(std::ostream& os, const basis_function& bf)
{
    os << "Eidx: " << bf.edge_index << ", ";
    os << "Midx: " << bf.matrix_index << ", ";
    os << "T-: " << bf.itminus << ", ";
    os << "T+: " << bf.itplus << ", " << bf.pminus << " " << bf.pplus;
    return os;
}

void make_function_space(const mesh&, std::vector<basis_function>&);

} // mommy
