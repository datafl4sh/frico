/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025-2026, Matteo Cicuttin - IV3IWE
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

#include "gmsh.h"
#include "constants.h"
#include "eigen.h"
#include "input_gmsh.h"
#include "geom_mesh.h"
#include "emw_solver.h"
#include "output_silo.h"
namespace frico::maxwell {

std::pair<ezvec3, ezvec3>
eval_fields(const simulation& sim, size_t ctx_number, const point& pt)
{
    const freq_context& context = sim.contexts[ctx_number];
    const mesh& msh = sim.msh;

    const double freq = context.frequency;
    const double omega = 2.0*M_PI*freq;
    const double k = omega*std::sqrt(MU0*EPS0);
    const std::complex<double> jomega{0.0, omega};
    const std::complex<double> jk{0.0, k};

    ezvec3 E = ezvec3::Zero();
    ezvec3 H = ezvec3::Zero();
    for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
        const auto& tri = msh.triangles[itri];
        const vec3 bar = barycenter(msh, tri);
        const edvec3 vR = (pt - bar).to_eigen();
        const double R = norm(pt - bar);
        const double _4piR = 4*M_PI*R;
        const double _4piR3 = _4piR*R*R;

        std::complex<double> jkR{0, k*R};
        std::complex<double> g = (std::exp(-jkR)/_4piR);
        ezvec3 J = context.tri_AJ.row(itri);
        std::complex<double> divJ = context.tri_AdivJ(itri);
        E += -jomega*MU0*(J - divJ*(1.0+jkR)*vR/(std::pow(k*R, 2)))*g;

        ezvec3 RcrossJ = vR.cross(J);
        std::complex<double> gradg = 
            -((1.0 + jkR)/_4piR3) * std::exp(jkR);
        H += RcrossJ * gradg;
    }

    return {E,H};
}

/**
 * @brief Compute fields E and H on all the points of the specified mesh
 * 
 * @param sim simulation
 * @param ctx_number context number
 * @param smpmsh Sampling mesh
 * @param E computed E-field
 * @param H computed H-field
 */
void eval_fields(const simulation& sim, size_t ctx_number,
    const mesh& smpmsh, zdfield& E, zdfield& H)
{
    E = zdfield::Zero(smpmsh.vertices.size(), 3);
    H = zdfield::Zero(smpmsh.vertices.size(), 3);

    #pragma omp parallel for
    for (size_t i = 0; i < smpmsh.vertices.size(); i++) {
        const auto& spt = smpmsh.vertices[i]; 
        auto [locE, locH] = eval_fields(sim, ctx_number, spt);
        E.row(i) = locE;
        H.row(i) = locH;
    }
}

/**
 * @brief Compute reflection coeff from impedance Z and system impedance Z0
 * 
 * @param Z Impedance
 * @param Z0 System impedance
 * @return std::complex<double> 
 */
std::complex<double>
compute_reflection_coefficient(std::complex<double> Z, double Z0)
{
    return (Z - Z0)/(Z + Z0);
}

/**
 * @brief Compute SWR from impedance Z and system impedance Z0
 * 
 * @param Z Impedance
 * @param Z0 System impedance
 * @return double 
 */
double
compute_swr(std::complex<double> Z, double Z0)
{
    std::complex<double> gamma = (Z - Z0)/(Z + Z0);
    return (1.0 + std::abs(gamma))/(1.0 - std::abs(gamma));
}

bool make_sampling_grid(frico::mesh& msh, const frico::point& c,
    double r, double h)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add("sampling");

    //int tagxy = gmsh::model::occ::addRectangle(c.x()-r, c.y()-r, c.z(), 2*r, 2*r);
    int tagyz = gmsh::model::occ::addRectangle(c.x()-r, c.y()-r, c.z(), 2*r, 2*r);
    //int tagxz = gmsh::model::occ::addRectangle(c.x()-r, c.y()-r, c.z(), 2*r, 2*r);
    gmsh::model::occ::rotate({{2, tagyz}},
        c.x(), c.y(), c.z(), c.x(), c.y()+1, c.z(), M_PI/2 );
    //gmsh::model::occ::rotate({{2, tagxz}},
    //    c.x(), c.y(), c.z(), c.x()+1, c.y(), c.z(), M_PI/2 );

    //gmsh::vectorpair out;
    //std::vector<gmsh::vectorpair> outmap;
    //gmsh::model::occ::fuse({ {2,tagxy}  }, {{2,tagyz}, {2,tagxz}}, out, outmap);

    gmsh::model::occ::synchronize();

    gmsh::vectorpair vp;
    gmsh::model::getEntities(vp);
    gmsh::model::mesh::setSize(vp, h);

    gmsh::model::mesh::generate(2);
    //gmsh::model::mesh::setOrder(1);

    frico::load_from_gmsh(msh, frico::load_mode::quick);

    gmsh::clear();
    gmsh::finalize();
    return true;
}

void
write_fields(const simulation& sim, size_t ctx_number)
{
    const auto& context = sim.contexts[ctx_number];

    std::string filename =
        sim.name + "_" + std::to_string(ctx_number) + ".silo";

    ddfield normals = ddfield::Zero(sim.msh.triangles.size(), 3);
    for (int i = 0; i < sim.msh.triangles.size(); i++) {
        normals.row(i) = normal(sim.msh, sim.msh.triangles[i]);
    }

    silo db;
    db.open(filename);
    db.add_mesh("mesh", sim.msh);
    db.add_variable("mesh", "J", context.tri_J, var_centering::zonal);
    db.add_variable("mesh", "normals", normals, var_centering::zonal);

    mesh smpmsh;
    make_sampling_grid(smpmsh, {0,0,0}, 1, 0.01);
    zdfield E = zdfield::Zero(smpmsh.vertices.size(), 3);
    zdfield H = zdfield::Zero(smpmsh.vertices.size(), 3);
    #pragma omp parallel for
    for (int i = 0; i < smpmsh.vertices.size(); i++) {
        const auto& pt = smpmsh.vertices[i];
        auto [locE, locH] = eval_fields(sim, ctx_number, pt);
        E.row(i) = locE;
        H.row(i) = locH;
    }

    db.add_mesh("smpmesh", smpmsh);
    db.add_variable("smpmesh", "E", E, var_centering::nodal);
    db.add_variable("smpmesh", "H", H, var_centering::nodal);
}

}