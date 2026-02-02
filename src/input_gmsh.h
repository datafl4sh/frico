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

#include <expected>
#include "geom_mesh.h"

namespace frico {

enum class meshing_error {
    gmsh_issue,
    lonely_edge,
    multiple_triangles,
};
struct meshing_error_info
{
    meshing_error   errtype;

    /* for bad_connectivity */
    int             tminus_tag;
    int             tplus_tag;
    int             offending_tag;
};

using merr_t = std::expected<bool, meshing_error_info>;

enum class load_mode {
    full,
    quick
};

/* Perhaps this went a bit too far... */
merr_t load_from_gmsh(mesh&);
merr_t load_from_gmsh(mesh&, const load_mode);
merr_t load_from_gmsh(mesh&, const load_mode, std::vector<int>&);
merr_t load_from_gmsh(const std::string&, mesh&);
merr_t load_from_gmsh(const std::string&, mesh&, const load_mode);
merr_t load_from_gmsh(const std::string&, mesh&, const std::vector<int>&);
merr_t load_from_gmsh(const std::string&, mesh&, const load_mode, const std::vector<int>&);

} // namespace frico