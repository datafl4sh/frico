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

#include "gmsh.h"
#include "input_gmsh.h"
#include "utils.h"
#include "sampling_planes.h"

namespace frico {

namespace priv {
enum class plane {
    invalid,
    xy,
    yz,
    xz
}; 

static plane whichplane(const std::string& str)
{
    if (str == "xy") return plane::xy;
    if (str == "yz") return plane::yz;
    if (str == "xz") return plane::xz;
    return plane::invalid;
}

} // namespace priv

static std::expected<sampling_plane, parse_error>
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

static std::expected<sampling_plane, parse_error>
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

std::expected<sampling_plane, parse_error>
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

bool make_sampling_plane_mesh(frico::mesh& msh, const sampling_plane& sp)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add("sampling");

    auto left = -sp.width/2.0;
    auto bottom = -sp.height/2.0;

    int tag = gmsh::model::occ::addRectangle(left, bottom, 0.0, sp.width, sp.height);

    gmsh::model::occ::rotate({{2, tag}},
        0.0, 0.0, 0.0, sp.normal.x(), sp.normal.y(), sp.normal.z(), M_PI/2 );

    gmsh::model::occ::synchronize();

    gmsh::vectorpair vp;
    gmsh::model::getEntities(vp);
    gmsh::model::mesh::setSize(vp, sp.h);

    gmsh::model::mesh::generate(2);

    frico::load_from_gmsh(msh, frico::load_mode::quick);

    gmsh::clear();
    gmsh::finalize();
    return true;
}

} // namespace frico