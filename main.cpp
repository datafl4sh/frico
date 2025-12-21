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

#include <iostream>
#include <cmath>

#include "gmsh.h"

#include "geom_mesh.h"
#include "input_gmsh.h"
#include "output_silo.h"
#include "quadratures.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Please specify GMSH .geo file" << std::endl;
        return 1;
    }

    try {
        gmsh::initialize(argc, argv);
        gmsh::open(argv[1]);
    }

    catch (const std::runtime_error& e) {
        std::cerr << "GMSH exception: " << e.what() << std::endl;
        return 1;
    }

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::setOrder(1);

    mommy::mesh msh;
    mommy::load_mesh_from_gmsh(msh);
    gmsh::finalize();

    std::cout << msh.vertices.size() << " " << msh.triangles.size() << std::endl;

    mommy::silo db("test.silo");
    db.add_mesh("mesh", msh);

    double l = 0.0;
    for (auto& be : msh.boundary_edges) {
        l += measure(msh, be);
    }
    std::cout << l << std::endl;

    double a = 0.0;
    for (auto& t : msh.triangles) {
        a += measure(msh, t);
    }
    std::cout << a << std::endl;

    std::vector<double> test;
    for (size_t i = 0; i < msh.triangles.size(); i++) {
        test.push_back(i);
    }

    db.add_variable("mesh", "test", test, mommy::var_centering::zonal);

    Eigen::VectorXd vals;
    vals = Eigen::VectorXd::Zero(msh.vertices.size());
    for (size_t i = 0; i < msh.vertices.size(); i++) {
        const auto& vtx = msh.vertices[i];
        auto val = std::sin(M_PI*vtx.x())*std::sin(M_PI*vtx.y());
        vals(i) = val;
    }
    db.add_variable("mesh", "vals", vals, mommy::var_centering::nodal);

    Eigen::Matrix<double, Eigen::Dynamic, 3> norms;
    norms = Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(msh.triangles.size(), 3);
    for (size_t i = 0; i < msh.triangles.size(); i++) {
        norms.row(i) = normal(msh, msh.triangles[i]);
    }
    db.add_variable("mesh", "normals", norms, mommy::var_centering::zonal);

    return 0;
}