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


#include "rwg_basis.h"

namespace frico {

void make_function_space(const mesh& msh, std::vector<basis_function>& bfs)
{
    bfs.reserve( num_internal_edges(msh) );
    size_t matrix_index = 0;
    for (size_t iedg = 0; iedg < msh.edge_neighbours.size(); iedg++) {
        const auto& en = msh.edge_neighbours[iedg];
        if (not en.itplus) {
            // If en.itplus is empty this
            // is a boundary edge: skip
            continue;
        }
    
        assert(en.itminus < msh.triangles.size());
        auto Tminus = msh.triangles[en.itminus];
        assert(en.itplus.value() < msh.triangles.size());
        auto Tplus = msh.triangles[en.itplus.value()];

        std::array<size_t, 3> ivtminus {Tminus.iv0, Tminus.iv1, Tminus.iv2};
        std::array<size_t, 3> ivtplus {Tplus.iv0, Tplus.iv1, Tplus.iv2};

        std::array<size_t, 3> evmap {2, 0, 1};
        assert( (en.loc_eminus < 3) and (en.loc_eplus.value() < 3) );
        size_t ipminus = ivtminus[evmap[en.loc_eminus]];
        size_t ipplus = ivtplus[evmap[en.loc_eplus.value()]];

        basis_function bf;
        bf.Aminus = measure(msh, Tminus);
        bf.Aplus = measure(msh, Tplus);
        bf.itminus = en.itminus;
        bf.itplus = en.itplus.value();
        assert(ipminus < msh.vertices.size());
        bf.pminus = msh.vertices[ipminus];
        assert(ipplus < msh.vertices.size());
        bf.pplus = msh.vertices[ipplus];
        assert(iedg < msh.edges.size());
        bf.length = measure(msh, msh.edges[iedg]);
        bf.edge_index = iedg;
        bf.matrix_index = matrix_index++;
        bf.interface = en.interface;
        bfs.push_back(bf);
    }
}

} // namespace frico