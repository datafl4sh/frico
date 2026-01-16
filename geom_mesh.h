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

#include <vector>
#include <array>
#include <optional>

#include <Eigen/Dense>

#include "geom_point.h"

namespace mommy {

struct edge {
    size_t  iv0;    // offset of first vertex
    size_t  iv1;    // offset of second vertex

    edge() = default;
    edge(size_t iv0p, size_t iv1p)
        : iv0(std::min(iv0p, iv1p)), iv1(std::max(iv0p, iv1p))
    {}

    bool operator<(const edge& other) const {
        assert(iv0 < iv1);
        assert(other.iv0 < other.iv1);
        bool a = iv0 < other.iv0;
        bool b = (iv0 == other.iv0) and (iv1 < other.iv1);
        return (a or b);
    }

    bool operator==(const edge& other) const {
        assert(iv0 < iv1);
        assert(other.iv0 < other.iv1);
        return (iv0 == other.iv0) and (iv1 == other.iv1);
    }
};

inline std::ostream&
operator<<(std::ostream& os, const edge& e)
{
    os << "Edg: (" << e.iv0 << ", " << e.iv1 << ") ";
    return os;
}

struct triangle {
    size_t  iv0;    // offset of first vertex
    size_t  iv1;    // offset of second vertex
    size_t  iv2;    // offset of third vertex
    int     tag;    // GMSH tag
};

inline std::ostream&
operator<<(std::ostream& os, const triangle& t)
{
    os << "Tri: (" << t.iv0 << ", " << t.iv1 << ", " << t.iv2 << "), " << t.tag;
    return os;
}

struct bedgeptr {
    size_t  offset;
    int     tag;
};

inline edge
deref(const std::vector<edge>& edges, bedgeptr bep)
{
    assert(bep.offset < edges.size());
    return edges[bep.offset];
}

struct neighbours {
    size_t                  itminus;    // T- index
    size_t                  loc_eminus; // Local edge number in T-
    std::optional<size_t>   itplus;     // T+ index
    std::optional<size_t>   loc_eplus;  // Local edge number in T+
    std::optional<size_t>   interface;  // Interface number, if interface
};

struct bf_info {
    size_t  edge_index;
    double  sign;
    point   p;
};

using triangle_bf_info = std::array<std::optional<bf_info>, 3>;

struct mesh {
    std::vector<point>          vertices;
    std::vector<bedgeptr>       boundary_edges;
    std::vector<edge>           edges;
    std::vector<triangle>       triangles;
    std::vector<neighbours>     edge_neighbours;
    std::vector<triangle_bf_info>   tbis;
};

template<typename T>
std::optional<size_t>
offset(const std::vector<T>& vec, const T& elem)
{
    auto itor = std::lower_bound(vec.begin(), vec.end(), elem);
    if (itor != vec.end() && *itor == elem) {
        return std::distance(vec.begin(), itor);
    }

    return {};
}

std::array<point, 2> points(const mesh&, const edge&);
std::array<point, 3> points(const mesh&, const triangle&);
point barycenter(const mesh&, const edge&);
point barycenter(const mesh&, const triangle&);
double measure(const mesh&, const edge&);
double measure(const mesh&, const triangle&);
bool is_boundary(const mesh&, const edge&);
size_t num_internal_edges(const mesh&);

inline std::array<edge, 3>
edges(const triangle& tri)
{
    return {{
        { tri.iv0, tri.iv1 },
        { tri.iv1, tri.iv2 },
        { tri.iv2, tri.iv0 }
    }};
}

Eigen::Matrix<double, 3, 1> normal(const mesh&, const triangle&);

}