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
std::complex<double> compute_reflection_coefficient(std::complex<double>, double);
double compute_swr(std::complex<double>, double);
bool make_sampling_grid(frico::mesh& msh, const frico::point& c,
    double r, double h);
void write_fields(const simulation& sim, size_t ctx_number);

namespace priv {
enum class plane {
    invalid,
    xy,
    yz,
    xz
};

inline plane whichplane(const std::string& str)
{
    if (str == "xy") return plane::xy;
    if (str == "yz") return plane::yz;
    if (str == "xz") return plane::xz;
    return plane::invalid;
}

}

inline std::expected<sampling_plane, parse_error>
parse_sampling_plane_short(const std::vector<std::string>& toks)
{
    assert(toks.size() == 4);

    sampling_plane ret;

    switch ( priv::whichplane(toks[0]) )
    {
        case priv::plane::xy:
            ret.normal = vec3{0.0, 0.0, 1.0};
            break;

        case priv::plane::yz:
            ret.normal = vec3{0.0, 1.0, 0.0};
            break;

        case priv::plane::xz:
            ret.normal = vec3{1.0, 0.0, 0.0};
            break;

        default:
            return std::unexpected(parse_error::out_of_range);
    }

    try {
        ret.width = std::stod(toks[1]);
        ret.height = std::stod(toks[2]);
        ret.h = std::stod(toks[3]);
    }
    catch (...) {
        return std::unexpected(parse_error::out_of_range);
    }

    return ret;
}

inline std::expected<sampling_plane, parse_error>
parse_sampling_plane_long(const std::vector<std::string>& toks)
{
    assert(toks.size() == 6);

    sampling_plane ret;

    try {
        double nx = std::stod(toks[0]);
        double ny = std::stod(toks[1]);
        double nz = std::stod(toks[2]);
        ret.normal = vec3{nx, ny, nz};
        ret.normal /= norm(ret.normal);
        ret.width = std::stod(toks[3]);
        ret.height = std::stod(toks[4]);
        ret.h = std::stod(toks[5]);
    }
    catch (...) {
        return std::unexpected(parse_error::out_of_range);
    }

    return ret;
}

inline std::expected<sampling_plane, parse_error>
parse_sampling_plane(const char *planespec)
{
    auto toks = split(planespec, ':');

    switch (toks.size()) {
        case 4:
            return parse_sampling_plane_short(toks);
            break;

        case 6:
            return parse_sampling_plane_long(toks);
            break;

        default:
            break;
    }

    return std::unexpected(parse_error::invalid_input);
}


bool make_sampling_plane_mesh(frico::mesh& msh, const sampling_plane& sp);

}

