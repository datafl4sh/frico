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

#include <complex>
#include <cassert>
#include <utility>
#include "output_silo.h"

namespace frico {

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
#ifdef SILO_USE_HDF5
    db_ = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
#else
    db_ = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
#endif
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
    if (not db_) {
        return false;
    }

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
    nodelist.reserve( 3*msh.triangles.size() + 2*msh.beams.size() );

    for (auto& t : msh.triangles) {
        nodelist.push_back( t.iv0 );
        nodelist.push_back( t.iv1 );
        nodelist.push_back( t.iv2 );
    }

    for (auto& b : msh.beams) {
        nodelist.push_back( b.iv0 );
        nodelist.push_back( b.iv1 );
    }

    int lnodelist = nodelist.size();
    int shapetype[] = { DB_ZONETYPE_TRIANGLE, DB_ZONETYPE_BEAM };
    int shapesize[] = {3, 2};
    int shapecounts[] = {
        static_cast<int>(msh.triangles.size()),
        static_cast<int>(msh.beams.size())
    };
    int nshapetypes = 2;
    int nnodes = msh.vertices.size();
    int nzones = msh.triangles.size() + msh.beams.size();
    int ndims = 3;

    std::string zlname = "zonelist_";
    zlname += name;

    if ( DBPutZonelist2(db_, zlname.c_str(), nzones, ndims,
        nodelist.data(), lnodelist, 0, 0, 0, shapetype, shapesize,
        shapecounts, nshapetypes, NULL) < 0 ) {
        return false;
    }

    if ( DBPutUcdmesh(db_, name.c_str(), ndims, NULL, coords, nnodes,
            nzones, zlname.c_str(), NULL, DB_DOUBLE, NULL) < 0 ) {
        return false;
    }

    return true;
}

bool
silo::add_variable(const std::string& mesh_name,
    const std::string& var_name,
    const std::vector<double>& var,
    var_centering centering)
{
    if (not db_) {
        return false;
    }

    switch (centering) {
        case var_centering::zonal: {
            if ( DBPutUcdvar1(db_, var_name.c_str(), mesh_name.c_str(),
                 var.data(), var.size(), NULL, 0, DB_DOUBLE,
                 DB_ZONECENT, NULL) < 0 ) {
                return false;
            }
        } break;

        case var_centering::nodal: {
            if ( DBPutUcdvar1(db_, var_name.c_str(), mesh_name.c_str(),
                 var.data(), var.size(), NULL, 0, DB_DOUBLE,
                 DB_NODECENT, NULL) < 0 ) {
                    return false;
            }
        } break;

        default:
            std::cerr << "Unreachable branch taken" << std::endl;
            std::terminate();
    }

    return true;
}

bool
silo::add_variable(const std::string& mesh_name,
    const std::string& var_name,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& var,
    var_centering centering)
{
    if (not db_) {
        return false;
    }

    switch (centering) {
        case var_centering::zonal: {
            if ( DBPutUcdvar1(db_, var_name.c_str(), mesh_name.c_str(),
                 var.data(), var.size(), NULL, 0, DB_DOUBLE,
                 DB_ZONECENT, NULL) < 0 ) {
                return false;
            }
        } break;

        case var_centering::nodal: {
            if ( DBPutUcdvar1(db_, var_name.c_str(), mesh_name.c_str(),
                 var.data(), var.size(), NULL, 0, DB_DOUBLE,
                 DB_NODECENT, NULL) < 0 ) {
                    return false;
            }
        } break;

        default:
            std::cerr << "Unreachable branch taken" << std::endl;
            std::terminate();
    }

    return true;
}

bool
silo::add_variable(const std::string& mesh_name,
    const std::string& var_name,
    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& var,
    var_centering centering)
{
    if (not db_) {
        return false;
    }

    const Eigen::Matrix<double, Eigen::Dynamic, 1> re = var.real();
    const Eigen::Matrix<double, Eigen::Dynamic, 1> im = var.imag();
    
    if ( not add_variable(mesh_name, var_name + "_re", re, centering) ) {
        return false;
    }
    if ( not add_variable(mesh_name, var_name + "_im", im, centering) ) {
        return false;
    }

    return true;
}


bool
silo::add_variable(const std::string& mesh_name,
    const std::string& var_name,
    const Eigen::Matrix<double, Eigen::Dynamic, 3>& var,
    var_centering centering)
{
    Eigen::VectorXd var_x = var.col(0);
    std::string vname_x = var_name + "_x";
    bool ok = add_variable(mesh_name, vname_x, var_x, centering);
    if (not ok) {
        return false;
    }

    Eigen::VectorXd var_y = var.col(1);
    std::string vname_y = var_name + "_y";
    ok = add_variable(mesh_name, vname_y, var_y, centering);
    if (not ok) {
        return false;
    }

    Eigen::VectorXd var_z = var.col(2);
    std::string vname_z = var_name + "_z";
    ok = add_variable(mesh_name, vname_z, var_z, centering);
    if (not ok) {
        return false;
    }

    std::string dbcurdir = curdir().value();
    std::string ename = dbcurdir + "/" + var_name;
    const char *names[] = { ename.c_str() };
    std::string def =  
        "{<" + dbcurdir + "/" + vname_x + ">," +
        "<" + dbcurdir + "/" + vname_y + ">," +
        "<" + dbcurdir + "/" + vname_z + ">}";
    const char *defs[] = { def.c_str() };
    int types[] = { DB_VARTYPE_VECTOR };
    std::string defname = dbcurdir + var_name + "defs";
    if ( DBPutDefvars(db_, defname.c_str(), 1, names, types, defs, NULL) < 0) {
        return false;
    }

    return true;
}

bool
silo::add_variable(const std::string& mesh_name,
    const std::string& var_name,
    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 3>& var,
    var_centering centering)
{
    /* Real part */
    Eigen::Matrix<double, Eigen::Dynamic, 3> real = var.real();
    std::string vname_re = var_name + "_real";
    if (not add_variable(mesh_name, vname_re, real, centering)) {
        return false;
    }
    
    /* Imaginary part */
    Eigen::Matrix<double, Eigen::Dynamic, 3> imag = var.imag();
    std::string vname_im = var_name + "_imag";
    if (not add_variable(mesh_name, vname_im, imag, centering)) {
        return false;
    }
    
    /* Magnitude */
    Eigen::Matrix<double, Eigen::Dynamic, 1> abs =
        Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(var.rows());
    for (Eigen::Index i = 0; i < var.rows(); i++) {
        Eigen::Matrix<std::complex<double>, 3, 1> r = var.row(i);
        abs(i) = std::sqrt( r.dot(r).real() );
    }
    std::string vname_mag = var_name + "_mag";
    if (not add_variable(mesh_name, vname_mag, abs, centering)) {
        return false;
    }

    return true;
}

bool
silo::add_curve(const std::string& name, const std::vector<double>& x,
    const std::vector<double>& y)
{
    if ( x.size() != y.size() ) {
        return false;
    }

    DBPutCurve(db_, name.c_str(), x.data(), y.data(),
        DB_DOUBLE, x.size(), nullptr);
    return true;
}

bool
silo::mkdir(const std::string& name)
{
    if (not db_) {
        return false;
    }
    DBMkDir(db_, name.c_str());
    return true;
}

bool
silo::chdir(const std::string& name)
{
    if (not db_) {
        return false;
    }
    DBSetDir(db_, name.c_str());
    return true;
}

std::optional<std::string>
silo::curdir(void) const
{
    if (not db_) {
        return {};
    }

    char dir[256];
    DBGetDir(db_, dir);
    return dir;
}

silo::~silo()
{
    close();
}

} // namespace frico