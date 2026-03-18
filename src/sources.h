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
#include <complex>

#include "eigen.h"
#include "geom_mesh.h"
#include "rwg_basis.h"

namespace frico {

/**************************************************************/
class plane_wave {
    using cd = std::complex<double>;

    edvec3  theta_hat;  // Theta unit vector (vertical)
    edvec3  phi_hat;    // Phi unit vector (horizontal)
    edvec3  dir;        // OPPOSITE propagation direction
    cd      Etheta;
    cd      Ephi;

public:
    plane_wave(double theta_i, double phi_i, cd Etheta_i, cd Ephi_i)
    {
        /* Compute the vector FROM the origin */
        dir = edvec3 {
            std::sin(theta_i) * std::cos(phi_i),
            std::sin(theta_i) * std::sin(phi_i),
            std::cos(theta_i)
        };

        /* Theta unit vector: rotating from Z towards XY plane */
        theta_hat = edvec3 {
            std::cos(theta_i) * std::cos(phi_i),
            std::cos(theta_i) * std::sin(phi_i),
            -std::sin(theta_i)
        };

        /* Phi unit vector: rotating from X towards YZ plane */
        phi_hat = edvec3 { -std::sin(phi_i), std::cos(phi_i), 0.0 };

        Etheta = Etheta_i;
        Ephi = Ephi_i;
    }

    ezvec3
    compute(double k, point pos) const
    {
        auto jk = std::complex<double>{0.0, k};
        ezvec3 E0 = Etheta*theta_hat + Ephi*phi_hat;
        /* Positive exponent: remember that dir is from (0,0,0) to pos */
        return E0 * std::exp(jk*dir.dot(pos.to_eigen()));  
    }
};

/**************************************************************/
struct delta_gap {
    std::string             name;
    int                     phys_entity;
    std::vector<int>        interfaces;
    std::complex<double>    voltage;
};


} // namespace frico