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
#include <fstream>

#include "gmsh.h"

#include "geom_mesh.h"
#include "input_gmsh.h"
#include "output_silo.h"
#include "quadratures.h"

namespace mommy {

struct basis_function
{
    double  Aminus;
    double  Aplus;
    size_t  itminus;
    size_t  itplus;
    point   pminus;
    point   pplus;
    std::vector<quadrature_point>   qminus;
    std::vector<quadrature_point>   qplus;
    double  length;
    size_t  edge_index;
    size_t  matrix_index;
};

std::ostream&
operator<<(std::ostream& os, const basis_function& bf)
{
    os << bf.edge_index << " " << bf.matrix_index << " " << bf.itminus;
    os << " " << bf.itplus;
    return os;
}

void populate_data(const mesh& msh, std::vector<basis_function>& bfs)
{
    bfs.reserve( num_internal_edges(msh) );
    size_t matrix_index = 0;
    for (size_t iedg = 0; iedg < msh.edge_neighbours.size(); iedg++) {
        const auto& en = msh.edge_neighbours[iedg];
        if (not en.itplus) {
            continue; // it is a boundary edge
        }

        size_t integration_degree = 4;

        auto Tminus = msh.triangles[en.itminus];
        auto Tplus  = msh.triangles[en.itplus.value()];

        std::array<size_t, 3> ivtminus {Tminus.iv0, Tminus.iv1, Tminus.iv2};
        std::array<size_t, 3> ivtplus {Tplus.iv0, Tplus.iv1, Tplus.iv2};

        std::array<size_t, 3> lvmap {2, 0, 1};
        size_t ipminus = ivtminus[lvmap[en.loc_eminus]];
        size_t ipplus = ivtplus[lvmap[en.loc_eplus.value()]];

        basis_function bf;
        bf.Aminus = measure(msh, Tminus);
        bf.Aplus = measure(msh, Tplus);
        bf.itminus = en.itminus;
        bf.itplus = en.itplus.value();
        bf.pminus = msh.vertices[ipminus];
        bf.pplus = msh.vertices[ipplus];
        bf.qminus = integrate(msh, Tminus, integration_degree);
        bf.qplus = integrate(msh, Tplus, integration_degree);
        bf.length = measure(msh, msh.edges[iedg]);
        bf.edge_index = iedg;
        bf.matrix_index = matrix_index++;
        bfs.push_back(bf);
    }
}

} //namespace mommy

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
        l += measure(msh, deref(msh.edges, be));
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

    std::cout << "Vertices: " << msh.vertices.size() << std::endl;
    std::cout << "Edges:    " << msh.edges.size() << std::endl;
    std::cout << "Cells:    " << msh.triangles.size() << std::endl;
    std::cout << "IntEdges: " << mommy::num_internal_edges(msh) << std::endl;

    std::vector<mommy::basis_function> bfs;
    mommy::populate_data(msh, bfs);

    for (size_t i = 0; i < bfs.size(); i++)
    {
        const auto& bf = bfs[i];
        std::string fname = "debug/basis_" + std::to_string(bf.edge_index) + "_minus.dat";
        std::ofstream ofs_minus(fname);
        for (auto& qpsminus : bf.qminus) {
            auto r = qpsminus.p;
            auto v = r - bf.pminus;
            auto rho = bf.length*v/(2.0*bf.Aminus);
            ofs_minus << r.x() << " " << r.y() << " " << v.x() << " " << v.y() << "\n";
        }

        fname = "debug/basis_" + std::to_string(bf.edge_index) + "_plus.dat";
        std::ofstream ofs_plus(fname);
        for (auto& qpsplus : bf.qplus) {
            auto r = qpsplus.p;
            auto v = bf.pplus - r;
            auto rho = bf.length*v/(2.0*bf.Aplus);
            ofs_plus << r.x() << " " << r.y() << " " << v.x() << " " << v.y() << "\n";
        }
    }

    return 0;
}

