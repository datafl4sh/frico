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

#include <print>
#include <cmath>

#include "eigen.h"


#include "geom_mesh.h"
#include "input_gmsh.h"
#include "output_silo.h"
#include "quadratures.h"
#include "rwg_basis.h"
#include "utils.h"
#include "sources.h"
#include "constants.h"
#include "emw_solver.h"
#include "emw_postpro_delta_gap.h"

int main(int argc, char **argv)
{
    //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    #ifndef _OPENMP
    std::print(
"FRICO v0.0 3D MoM solver - IV3IWE Matteo Cicuttin (C) 2025-2026\n\n"
    );
    #else
    std::print(
"FRICO v0.0 3D MoM solver - IV3IWE Matteo Cicuttin (C) 2025-2026 [OpenMP]\n\n"
    );
    #endif

    frico::maxwell::simulation sim;

    const char *arg_source = nullptr;
    const char *arg_geo_path = nullptr;
    const char *arg_skiptags = nullptr;
    const char *arg_frequency = nullptr;
    const char *arg_range_expr = nullptr;
    const char *arg_simname = "default";

    int opt;
    while ((opt = getopt(argc, argv, "Adf:g:k:n:s:SR:x:Z:")) != -1) {
        switch (opt) {
        case 'A':
            sim.cfg.approx_matrix = true;
            break;
        case 'd':
            sim.cfg.dump_matrices = true;
            break;
        case 'f':
            arg_frequency = optarg;
            break;
        case 'g':
            arg_geo_path = optarg;
            break;
        case 'k':
            sim.cfg.degree = std::stoull(optarg);
            break;
        case 'n':
            arg_simname = optarg;
            break;
        case 's':
            arg_source = optarg;
            break;
        case 'S':
            sim.cfg.force_symmetry = true;
            break;
        case 'R':
            arg_range_expr = optarg;
            break;
        case 'x':
            arg_skiptags = optarg;
            break;
        case 'Z':
            sim.cfg.Z0 = std::stod(optarg);
            break;

        default:
            std::println(stderr, "Invalid argument");
            return EXIT_FAILURE;
        }
    }

    /* (1) Check if geometry was specified */
    if (arg_geo_path == nullptr) {
        std::println(stderr, "No geometry specified (-g)");
        return EXIT_FAILURE;
    }

    /* (2) Check if frequency was specified, either single or sweep */
    auto opt_freqs = frico::parse_frequency_parameters(arg_frequency, arg_range_expr);
    if (not opt_freqs) {
        return EXIT_FAILURE;
    }

    /* (3) If there is a list of surfaces to skip, process it */
    if (arg_skiptags) {
        auto exp_skiptags = frico::parse_integer_list(arg_skiptags);
        if (exp_skiptags.has_value()) {
            sim.skiptags = *exp_skiptags;
        } else {
            std::println(stderr, "Error parsing the argument of -x");
            return EXIT_FAILURE;
        }
    }

    if ( not frico::maxwell::init_simulation(sim, arg_simname, arg_geo_path) ) {
        return EXIT_FAILURE;
    }

    if (not arg_source) {
        std::println(stderr, "No sources specified. Exiting.");
        return EXIT_FAILURE;
    }

    const auto& pgs = sim.msh.physgroups;
    if ( pgs.find(arg_source) == pgs.end() ) {
        std::println(stderr, 
            "Unknown physical group \"{}\": cannot enable source",
            arg_source);
        return EXIT_FAILURE;
    }
    const auto& pg = sim.msh.physgroups[arg_source];
    frico::delta_gap dg;
    dg.name = pg.name;
    dg.phys_entity = pg.tag;
    dg.interfaces = pg.entityTags;
    dg.voltage = 1.0;


    frico::maxwell::init_sweep(sim, *opt_freqs);
    frico::maxwell::do_sweep(sim, dg);

    return EXIT_SUCCESS;
}
