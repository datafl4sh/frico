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

#include "output_silo.h"

namespace mommy {

silo::silo()
    : db_(nullptr)
{}

silo::silo(const std::string& filename)
{
    if (not open(filename)) {
        throw std::runtime_error("Can't open " + filename);
    }
}

silo::silo(const char *filename)
{
    if (not open(filename)) {
        throw std::runtime_error("Can't open " + std::string(filename));
    }
}

bool
silo::open(const std::string& filename)
{
    return open(filename.c_str());
}

bool
silo::open(const char *filename)
{
    db_ = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
    if (db_ == nullptr) {
        std::cerr << "SILO: Can't open " << filename << "\n";
        return false;
    }

    return true;
}

bool
silo::is_open() const
{
    return db_ != nullptr;
}

void
silo::close()
{
    if (db_) {
        DBClose(db_);
    }

    db_ = nullptr;
}

bool
silo::add_mesh(const std::string& name, const mesh& msh)
{
    std::vector<double> x_coords, y_coords, z_coords;
    x_coords.reserve(msh.vertices.size());
    y_coords.reserve(msh.vertices.size());
    z_coords.reserve(msh.vertices.size());

    for (auto& vtx : msh.vertices) {
        x_coords.push_back( vtx.x() );
        y_coords.push_back( vtx.y() );
        z_coords.push_back( vtx.z() );
    }

    double *coords[] = {x_coords.data(), y_coords.data(), z_coords.data() };

    std::vector<int> nodelist;
    nodelist.reserve( 3*msh.triangles.size() );

    for (auto& t : msh.triangles)
    {
        nodelist.push_back( t.iv0 + 1 );
        nodelist.push_back( t.iv1 + 1 );
        nodelist.push_back( t.iv2 + 1 );
    }

    int lnodelist = nodelist.size();

    int shapetype[] = { DB_ZONETYPE_TRIANGLE };
    int shapesize[] = {3};
    int shapecounts[] = { static_cast<int>(msh.triangles.size()) };
    int nshapetypes = 1;
    int nnodes = msh.vertices.size();
    int nzones = msh.triangles.size();
    int ndims = 3;

    std::string zlname = "zonelist_";
    zlname += name;

    if ( DBPutZonelist2(db_, zlname.c_str(), nzones, ndims,
        nodelist.data(), lnodelist, 1, 0, 0, shapetype, shapesize,
        shapecounts, nshapetypes, NULL) < 0 ) {
        std::cerr << "DBPutZoneList2() failed\n";
        return false;
    }

    std::cout << name.c_str() << std::endl;

    if ( DBPutUcdmesh(db_, name.c_str(), ndims, NULL, coords, nnodes,
            nzones, zlname.c_str(), NULL, DB_DOUBLE, NULL) < 0 ) {
        std::cerr << "DBPutUcdmesh() failed\n";
        return false;
    }

    return true;
}

silo::~silo()
{
    close();
}

} // namespace mommy