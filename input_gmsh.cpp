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

#include <vector>
#include <algorithm>
#include <optional>
#include <cassert>

#include "geom_mesh.h"
#include "gmsh.h"

namespace mommy {

static void
gmsh_get_vertices(mesh& msh, std::vector<std::optional<size_t>>& node_tag2ofs)
{
    std::vector<size_t>     nodeTags;
    std::vector<double>     coords;
    std::vector<double>     paraCoords;

    gmsh::model::mesh::getNodes(nodeTags, coords, paraCoords, -1, -1, true, false);

    auto maxtag_pos = std::max_element(nodeTags.begin(), nodeTags.end());
    
    size_t maxtag = 0;
    if (maxtag_pos != nodeTags.end())
        maxtag = (*maxtag_pos)+1;

    msh.vertices.resize( 3*nodeTags.size() );

    for (size_t i = 0; i < nodeTags.size(); i++) {
        auto x = coords[3*i+0];
        auto y = coords[3*i+1];
        auto z = coords[3*i+2];
        msh.vertices[i] = point(x,y,z);
    }

    node_tag2ofs.resize(maxtag);
    for (size_t i = 0; i < nodeTags.size(); i++) {
        node_tag2ofs.at( nodeTags[i] ) = i;
    }
}

static void
gmsh_get_triangles(mesh& msh, const std::vector<std::optional<size_t>>& node_tag2ofs)
{
    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities, 2/*dimension*/);
    size_t subdom_id = 0;
    for (auto [dim, tag] : entities)
    {
        assert(dim == 2);
        std::vector<int> elemTypes;
        gmsh::model::mesh::getElementTypes(elemTypes, dim, tag);
        for (auto& elemType : elemTypes)
        {
            std::vector<size_t> elemTags;
            std::vector<size_t> elemNodeTags;
            gmsh::model::mesh::getElementsByType(elemType, elemTags, elemNodeTags, tag);
            auto nodesPerElem = elemNodeTags.size()/elemTags.size();
            assert( elemTags.size() * nodesPerElem == elemNodeTags.size() );
        
            msh.triangles.reserve( 3*elemTags.size() );

            for (size_t i = 0; i < elemTags.size(); i++)
            {
                auto base = nodesPerElem * i;

                auto node0_tag = elemNodeTags[base + 0];
                assert(node0_tag < node_tag2ofs.size());
                auto node0_ofs = node_tag2ofs[node0_tag].value(); //must be valid

                auto node1_tag = elemNodeTags[base + 1];
                assert(node1_tag < node_tag2ofs.size());
                auto node1_ofs = node_tag2ofs[node1_tag].value(); //must be valid

                auto node2_tag = elemNodeTags[base + 2];
                assert(node2_tag < node_tag2ofs.size());
                auto node2_ofs = node_tag2ofs[node2_tag].value(); //must be valid

                triangle t{node0_ofs, node1_ofs, node2_ofs, tag};
                msh.triangles.push_back(t);

                msh.edges.push_back( {node0_ofs, node1_ofs} );
                msh.edges.push_back( {node1_ofs, node2_ofs} );
                msh.edges.push_back( {node2_ofs, node0_ofs} );
            }
        }

        subdom_id++;
    }

    std::sort(msh.edges.begin(), msh.edges.end());
    msh.edges.erase(
        std::unique(msh.edges.begin(), msh.edges.end()),
        msh.edges.end()
    );
}

static void
gmsh_get_boundary_edges(mesh& msh, const std::vector<std::optional<size_t>>& node_tag2ofs)
{
    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities, 1/*dimension*/);
    size_t subdom_id = 0;
    for (auto [dim, tag] : entities)
    {
        assert(dim == 1);
        std::vector<int> elemTypes;
        gmsh::model::mesh::getElementTypes(elemTypes, dim, tag);
        for (auto& elemType : elemTypes)
        {
            std::vector<size_t> elemTags;
            std::vector<size_t> elemNodeTags;
            gmsh::model::mesh::getElementsByType(elemType, elemTags, elemNodeTags, tag);
            auto nodesPerElem = elemNodeTags.size()/elemTags.size();
            assert( elemTags.size() * nodesPerElem == elemNodeTags.size() );
        
            msh.boundary_edges.reserve( 2*elemTags.size() );

            for (size_t i = 0; i < elemTags.size(); i++)
            {
                auto base = nodesPerElem * i;

                auto node0_tag = elemNodeTags[base + 0];
                assert(node0_tag < node_tag2ofs.size());
                auto node0_ofs = node_tag2ofs[node0_tag].value(); //must be valid

                auto node1_tag = elemNodeTags[base + 1];
                assert(node1_tag < node_tag2ofs.size());
                auto node1_ofs = node_tag2ofs[node1_tag].value(); //must be valid

                edge e{node0_ofs, node1_ofs};

                bedgeptr bep;
                bep.offset = offset(msh.edges, e).value(); // must exist
                bep.tag = tag;
                msh.boundary_edges.push_back(bep);
            }
        }

        subdom_id++;
    }
}

void
load_mesh_from_gmsh(mesh& msh)
{
    std::vector<std::optional<size_t>> node_tag2ofs;
    gmsh_get_vertices(msh, node_tag2ofs);
    gmsh_get_triangles(msh, node_tag2ofs);
    gmsh_get_boundary_edges(msh, node_tag2ofs);

    size_t max_index = 0;
    for (auto& t : msh.triangles) {
        max_index = std::max( {max_index, t.iv0, t.iv1, t.iv2} );
    }

    std::vector<std::optional<size_t>> used(max_index+1);
    for (auto& t : msh.triangles) {
        used[t.iv0] = 1;
        used[t.iv1] = 1;
        used[t.iv2] = 1;
    }

    int ci = 0;
    for (auto& uopt : used) {
        if (uopt) {
            uopt = ci++;
        }
    }

    std::vector<point> vertices(ci);
    for (size_t i = 0; i < used.size(); i++) {
        const auto& uopt = used[i];
        if (uopt) {
            vertices[uopt.value()] = msh.vertices[i];
        }
    }
    std::swap(vertices, msh.vertices);
    vertices.clear();

    for (auto& t : msh.triangles) {
        t.iv0 = used[t.iv0].value();
        t.iv1 = used[t.iv1].value();
        t.iv2 = used[t.iv2].value();
    }

    for (auto& e : msh.edges) {
        e.iv0 = used[e.iv0].value();
        e.iv1 = used[e.iv1].value();
    }
}

} // namespace mommy