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

#include <memory>

#include "geom_mesh.h"
#include "rwg_basis.h"
#include "sources.h"
#include "utils.h"

namespace frico::maxwell {

struct config
{
    size_t          degree = 1;
    bool            approx_matrix = false;
    bool            force_symmetry = false;
    bool            dump_matrices = false;
};
struct freq_context {
    size_t                      ctx_number; // Incremental context number
    double                      frequency;  // Frequency
    zdmatrix                    Z;          // Impedance matrix
    zdvector                    V;          // Voltage vector (rhs)
    zdvector                    I;          // Current vector (unknown)
    zdfield                     tri_J;      // Tri-by-tri current density
    zdfield                     tri_AJ;     
    zdvector                    tri_AdivJ;
    std::complex<double>        fp_P;
    std::complex<double>        fp_Z;
    double                      gain;
};

enum class simulation_type {
    antenna,
    radar
};

struct simulation {
    std::string                 name;       // Simulation name
    simulation_type             type;       // Type of simulation
    config                      cfg;        // Simulation config
    mesh                        msh;        // Mesh
    std::vector<int>            skiptags;   // Surfaces to skip
    std::vector<basis_function> bfuncs;     // Basis functions
    std::vector<freq_context>   contexts;   // Data for each frequency
    std::unique_ptr<excitation> excit;
};

bool init_simulation(simulation&, const std::string&,
    const std::string&);
bool run(simulation&);
bool init_sweep(simulation&, const frequency_range&);

} //namespace frico::maxwell
