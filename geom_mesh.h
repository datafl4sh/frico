#pragma once

#include <vector>

#include "geom_point.h"

namespace mommy {

struct edge {
    size_t  iv0;    // offset of first vertex
    size_t  iv1;    // offset of second vertex
    int     tag;    // GMSH tag
};

inline std::ostream&
operator<<(std::ostream& os, const edge& e)
{
    os << "Edg: (" << e.iv0 << ", " << e.iv1 << "), " << e.tag;
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

struct mesh {
    std::vector<point>      vertices;
    std::vector<edge>       boundary_edges;
    std::vector<edge>       edges;
    std::vector<triangle>   triangles;
};

}