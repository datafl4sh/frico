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

#include <complex>

#include "eigen.h"
#include "geom_mesh.h"
#include "rwg_basis.h"

namespace frico {

/**************************************************************/
class excitation {
public:
    excitation(){}
    virtual std::complex<double>
    compute(const mesh&, const basis_function &) const = 0;
    virtual ~excitation(){};
};

/**************************************************************/
class plane_wave : public excitation {
    edvec3 E0_;
    edvec3 kinc_;

public:
    plane_wave(const edvec3& E0, const edvec3 kinc)
        : E0_(E0), kinc_(kinc) {};
    std::complex<double>
    compute(const mesh&, const basis_function &) const final override;
};

/**************************************************************/
class delta_gap : public excitation {
    int                     interface_;
    std::complex<double>    voltage_;
    double                  omega_;

public:
    delta_gap(int interface, std::complex<double>& voltage, double omega)
        : interface_(interface), voltage_(voltage), omega_(omega)
    {}

    std::complex<double>
    compute(const mesh& msh, const basis_function &bf) const final override;
};


} // namespace frico