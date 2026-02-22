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

#include <print>
#include <chrono>

#include "eigen.h"
#include <highfive/H5Easy.hpp>

#include "geom_mesh.h"
#include "rwg_basis.h"
#include "sources.h"
#include "utils.h"
#include "lapack_wrappers.h"

namespace frico::maxwell {

struct config
{
    size_t          degree = 1;
    bool            approx_matrix = false;
    bool            force_symmetry = false;
    bool            dump_matrices = false;
    bool            force_reorient_deltagap = false;
    bool            verbose = false;
    double          Z0 = 50.0;
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
    std::vector<double>         delta_gap_signs;
};

bool init_simulation(simulation&, const std::string&,
    const std::string&);

void compute_matrix(simulation&, size_t);
void compute_matrix_approx(simulation&, size_t);
void update_rhs(simulation&, size_t, const delta_gap&);
void update_rhs(simulation&, size_t, const plane_wave&);

bool reorient_deltagap_edges(const simulation&, const delta_gap&, std::vector<double>&);

template<typename Source>
bool run_context(simulation& sim, size_t ctx_number, const Source& src)
{
    freq_context& context = sim.contexts[ctx_number];

    std::println("********************************************");
    std::println("Sweep step {}: {} Hz", ctx_number, context.frequency);

    auto system_size = num_internal_edges(sim.msh);
    context.Z = zdmatrix::Zero(system_size, system_size);
    context.V = zdvector::Zero(system_size);

    std::print("  Assemblying linear system..."); std::fflush(stdout);
    const auto asm_start{std::chrono::steady_clock::now()};
    if (sim.cfg.approx_matrix) {
        compute_matrix_approx(sim, ctx_number);
    } else {
        compute_matrix(sim, ctx_number);
    }

    update_rhs(sim, ctx_number, src);

    const auto asm_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
    std::println("{} seconds", asm_elapsed_seconds);

    if (sim.cfg.force_symmetry) {
        context.Z = ((context.Z+context.Z.transpose())/2.0).eval();
    }

    std::print("  Solving linear system..."); std::fflush(stdout);
    const auto start{std::chrono::steady_clock::now()};
    context.I = context.Z.lu().solve(context.V);
    //solve_general(context.Z, context.V, context.I);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::println("{} seconds", elapsed_seconds);

    context.tri_AJ = zdfield::Zero(sim.msh.triangles.size(), 3);
    context.tri_AdivJ = zdvector::Zero(sim.msh.triangles.size());
    context.tri_J = zdfield::Zero(sim.msh.triangles.size(), 3);

    for (const auto& bf : sim.bfuncs) {
        const auto& Tminus = sim.msh.triangles[bf.itminus];
        const auto& Tplus = sim.msh.triangles[bf.itplus];

        auto bar_Tminus = barycenter(sim.msh, Tminus);
        auto bar_Tplus = barycenter(sim.msh, Tplus);

        auto A_Tminus = measure(sim.msh, Tminus);
        auto A_Tplus = measure(sim.msh, Tplus);

        std::complex<double> Iedge = context.I(bf.matrix_index);

        ezvec3 Jminus = bf.eval_minus(bar_Tminus) * Iedge;
        ezvec3 Jplus = bf.eval_plus(bar_Tplus) * Iedge;

        context.tri_AJ.row(bf.itminus) += A_Tminus * Jminus;
        context.tri_AJ.row(bf.itplus) += A_Tplus * Jplus;

        context.tri_AdivJ(bf.itminus) += A_Tminus * bf.div_minus(bar_Tminus) * Iedge;
        context.tri_AdivJ(bf.itplus) += A_Tplus * bf.div_plus(bar_Tplus) * Iedge;

        context.tri_J.row(bf.itminus) += Jminus;
        context.tri_J.row(bf.itplus) += Jplus;
    }

    return true;
}





void postpro_context(simulation&, size_t, const plane_wave&);
void write_file_headers(const simulation&, const plane_wave&);


bool init_sweep(simulation&, const frequency_range&);


template<typename Source>
bool do_sweep(simulation& sim, const Source& src)
{
    write_file_headers(sim, src);
    reorient_deltagap_edges(sim, src, sim.delta_gap_signs);

    for (size_t ctx_num = 0; ctx_num < sim.contexts.size(); ctx_num++) {
        run_context(sim, ctx_num, src);
        postpro_context(sim, ctx_num, src);

        /* Dump matrices, if needed*/
        if (sim.cfg.dump_matrices) {
            const auto& context = sim.contexts[ctx_num];
            std::string h5fn = "frico_" + std::to_string(ctx_num) + ".h5";
            H5Easy::File file(h5fn, H5Easy::File::Truncate);
            file.createDataSet("/frico/Z", context.Z);
            file.createDataSet("/frico/V", context.V);
        }

        /* dealloc matrix when we're done */
        sim.contexts[ctx_num].Z.resize(0,0);
    }
    return true;
}



} //namespace frico::maxwell
