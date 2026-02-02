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

#include <vector>
#include <algorithm>
#include <optional>
#include <cassert>
#include <expected>

#include "geom_mesh.h"
#include "input_gmsh.h"
#include "gmsh.h"

namespace frico {

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

static inline bool
skip_tag(const std::vector<int>& tags, int qtag)
{
    if (tags.size() > 32) { /* YES, I benchmarked it. */
        return std::binary_search(tags.begin(), tags.end(), qtag);
    }

    for (const auto& tag : tags) {
        if (tag == qtag) {
            return true;
        }
    }
    return false;
}

static void
gmsh_get_triangles(mesh& msh, const std::vector<std::optional<size_t>>& node_tag2ofs,
    const std::vector<int>& skiptags)
{
    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities, 2/*dimension*/);

    for (auto [dim, tag] : entities)
    {
        if ( skip_tag(skiptags, tag) ) {
            continue;
        }

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
                assert(node_tag2ofs[node0_tag]);
                auto node0_ofs = *node_tag2ofs[node0_tag];

                auto node1_tag = elemNodeTags[base + 1];
                assert(node1_tag < node_tag2ofs.size());
                assert(node_tag2ofs[node1_tag]);
                auto node1_ofs = *node_tag2ofs[node1_tag];

                auto node2_tag = elemNodeTags[base + 2];
                assert(node2_tag < node_tag2ofs.size());
                assert(node_tag2ofs[node2_tag]);
                auto node2_ofs = *node_tag2ofs[node2_tag];

                triangle t{node0_ofs, node1_ofs, node2_ofs, tag};
                msh.triangles.push_back(t);

                msh.edges.push_back( {node0_ofs, node1_ofs} );
                msh.edges.push_back( {node1_ofs, node2_ofs} );
                msh.edges.push_back( {node2_ofs, node0_ofs} );
            }
        }
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
                assert(node_tag2ofs[node0_tag]);
                auto node0_ofs = *node_tag2ofs[node0_tag];

                auto node1_tag = elemNodeTags[base + 1];
                assert(node1_tag < node_tag2ofs.size());
                assert(node_tag2ofs[node1_tag]);
                auto node1_ofs = *node_tag2ofs[node1_tag];

                edge e{node0_ofs, node1_ofs};

                auto opt_ofs = offset(msh.edges, e);
                if (not opt_ofs) {
                    std::cout << "Edge not found: " << e << " " << tag << std::endl;
                }
                assert(opt_ofs);
                bedgeptr bep;
                bep.offset = *opt_ofs;
                bep.tag = tag;
                msh.boundary_edges.push_back(bep);
            }
        }
    }
}

static std::expected<bool, meshing_error_info>
compute_connectivity(mesh& msh)
{
    std::vector<int> flags(msh.edges.size(), 0);

    msh.edge_neighbours.resize( msh.edges.size() );

    for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
        const auto& tri = msh.triangles[itri];
        auto edgs = edges(tri);
        std::array<size_t, 3> ofss;
        ofss[0] = offset(msh.edges, edgs[0]).value();
        ofss[1] = offset(msh.edges, edgs[1]).value();
        ofss[2] = offset(msh.edges, edgs[2]).value();

        for (int iedg = 0; iedg < 3; iedg++) {
            auto ofs = ofss[iedg];
            auto& en = msh.edge_neighbours[ofs];
            if (en.itplus) {
                meshing_error_info mei;
                mei.errtype = meshing_error::bad_connectivity;
                mei.tminus_tag = msh.triangles[en.itminus].tag;
                mei.tplus_tag = msh.triangles[*en.itplus].tag;
                mei.offending_tag = msh.triangles[itri].tag;
                return std::unexpected(mei);
            }
            assert(not en.loc_eplus);
            if (flags[ofs]) {
                if ( en.itminus < itri) {
                    en.itplus = itri;
                    en.loc_eplus = iedg;
                } else {
                    en.itplus = en.itminus;
                    en.itminus = itri;
                    en.loc_eplus = en.loc_eminus;
                    en.loc_eminus = iedg;
                }
            } else {
                en.itminus = itri;
                en.loc_eminus = iedg;
                flags[ofs] = 1;
            }
        
            // valid(T+) => ( T- < T+ )
            assert( (not en.itplus) or (en.itminus < en.itplus.value()) );
        }
    }

    for (auto& be : msh.boundary_edges) {
        assert(be.offset < msh.edge_neighbours.size());
        if (msh.edge_neighbours[be.offset].itplus) {
            msh.edge_neighbours[be.offset].interface = be.tag;
        }
    }

    /* to obtain basis functions from triangles */
    for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
        const auto& tri = msh.triangles[itri];
        auto edgs = edges(tri);
        std::array<size_t, 3> ofss;
        ofss[0] = offset(msh.edges, edgs[0]).value();
        ofss[1] = offset(msh.edges, edgs[1]).value();
        ofss[2] = offset(msh.edges, edgs[2]).value();

        triangle_bf_info tbi;

        std::array<size_t, 3> evmap {2, 0, 1};

        for (int iedg = 0; iedg < 3; iedg++) {
            auto ofs = ofss[iedg];
            auto& en = msh.edge_neighbours[ofs];
            if (not en.itplus) {
                continue;
            }
        
            /* negative */
            if (itri == en.itminus) {
                auto T = msh.triangles[en.itminus];
                std::array<size_t, 3> ivt{T.iv0, T.iv1, T.iv2};
                auto ip = ivt[evmap[iedg]];
                
                bf_info bi;
                bi.edge_index = ofss[iedg];
                bi.sign = -1;
                bi.p = msh.vertices[ip];
                tbi[iedg] = bi;
            }

            /* positive */
            if (itri == en.itplus.value()) {
                auto T = msh.triangles[en.itplus.value()];
                std::array<size_t, 3> ivt{T.iv0, T.iv1, T.iv2};
                auto ip = ivt[evmap[iedg]];
                
                bf_info bi;
                bi.edge_index = ofss[iedg];
                bi.sign = +1;
                bi.p = msh.vertices[ip];
                tbi[iedg] = bi;
            }
        }

        msh.tbis.push_back(tbi);
    }

    return true;
}

merr_t
load_from_gmsh(mesh& msh, const load_mode mode, const std::vector<int>& skiptags)
{
    std::vector<std::optional<size_t>> node_tag2ofs;
    gmsh_get_vertices(msh, node_tag2ofs);
    gmsh_get_triangles(msh, node_tag2ofs, skiptags);
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
            vertices[*uopt] = msh.vertices[i];
        }
    }
    std::swap(vertices, msh.vertices);
    vertices.clear();

    for (auto& t : msh.triangles) {
        assert(used[t.iv0]);
        t.iv0 = *used[t.iv0];
        assert(used[t.iv1]);
        t.iv1 = *used[t.iv1];
        assert(used[t.iv2]);
        t.iv2 = *used[t.iv2];
    }

    for (auto& e : msh.edges) {
        assert(used[e.iv0]);
        e.iv0 = *used[e.iv0];
        assert(used[e.iv1]);
        e.iv1 = *used[e.iv1];
    }

    if (mode == load_mode::full) {
        auto ccret = compute_connectivity(msh);
        if (not ccret) {
            return ccret;
        }
    }

    return true;
}

merr_t
load_from_gmsh(mesh& msh)
{
    return load_from_gmsh(msh, load_mode::full, {});
}

merr_t
load_from_gmsh(mesh& msh, const load_mode mode)
{
    return load_from_gmsh(msh, mode, {});
}

merr_t
load_from_gmsh(const std::string& filename, mesh& msh)
{
    return load_from_gmsh(filename, msh, load_mode::full, {});
}

merr_t
load_from_gmsh(const std::string& filename, mesh& msh, const load_mode mode)
{
    return load_from_gmsh(filename, msh, mode, {});
}

merr_t
load_from_gmsh(const std::string& filename,
    mesh& msh, const std::vector<int>& skiptags)
{
    return load_from_gmsh(filename, msh, load_mode::full, skiptags);
}

merr_t
load_from_gmsh(const std::string& filename, mesh& msh,
    const load_mode mode, const std::vector<int>& skiptags)
{
    try {
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", 1);
        gmsh::open(filename);
    }

    catch (const std::runtime_error& e) {
        std::cerr << "GMSH exception: " << e.what() << std::endl;
        return false;
    }

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::setOrder(1);

    merr_t ret = load_from_gmsh(msh, mode, skiptags);
    gmsh::clear();
    gmsh::finalize();

    return ret;
}


} // namespace frico