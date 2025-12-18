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

#include <cmath>

#include "geom_mesh.h"

namespace mommy {

std::array<point, 2>
points(const mesh& msh, const edge& e)
{
    return {
        msh.vertices[e.iv0],
        msh.vertices[e.iv1]
    };
}

point
barycenter(const mesh& msh, const edge& e)
{
    auto pts = points(msh, e);
    return (pts[0] + pts[1])/2.0;
}

double
measure(const mesh& msh, const edge& e)
{
    auto pts = points(msh, e);
    auto d = pts[1] - pts[0];
    return std::hypot(d.x(), d.y());
}

std::array<point, 3>
points(const mesh& msh, const triangle& t)
{
    return {
        msh.vertices[t.iv0],
        msh.vertices[t.iv1],
        msh.vertices[t.iv2]
    };
}

point
barycenter(const mesh& msh, const triangle& t)
{
    auto pts = points(msh, t);
    return (pts[0] + pts[1] + pts[2])/3.0;
}

double
measure(const mesh& msh, const triangle& t)
{
    auto pts = points(msh, t);
    auto e0 = pts[1] - pts[0];
    auto e1 = pts[2] - pts[1];
    return 0.5*norm(cross(e0,e1));
}

} // namespace mommy