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
#include <chrono>

#include "eigen.hpp"

#include <highfive/H5Easy.hpp>

#include "gmsh.h"

#include "geom_mesh.h"
#include "input_gmsh.h"
#include "output_silo.h"
#include "quadratures.h"
#include "rwg_basis.hpp"

#define MU0     1.256637061435917e-06
#define EPS0    8.8541878188e-12
namespace mommy {

struct config
{
    double  frequency;
    size_t  degree;
};

ddvector prjtest(const mesh& msh, const std::vector<basis_function>& bfs)
{
    ddmatrix mass = ddmatrix::Zero(bfs.size(), bfs.size());
    ddvector rhs = ddvector::Zero(bfs.size());

    std::vector<size_t> cmap;
    cmap.resize( msh.edges.size() );
    for (const auto& bf : bfs) {
        cmap[bf.edge_index] = bf.matrix_index;
    }

    for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
        const auto& tri = msh.triangles[itri];
        const auto& tbi = msh.tbis[itri];

        auto qps = integrate(msh, tri, 2);
        for (const auto& qp : qps) {
            for (size_t i = 0; i < 3; i++) {
                if (not tbi[i]) { continue; }
                auto bi = *tbi[i];
                auto gi = cmap[ bi.edge_index ];
                const auto& ibf = bfs[ gi ];
                edvec3 ibfval = (bi.sign > 0) ? ibf.eval_plus(qp.p) : ibf.eval_minus(qp.p);

                for (size_t j = 0; j < 3; j++) {
                    if (not tbi[j]) { continue; }
                    auto bj = *tbi[j];
                    auto gj = cmap[ bj.edge_index ];
                    const auto& jbf = bfs[ gj ];
                    edvec3 jbfval = (bj.sign > 0) ? jbf.eval_plus(qp.p) : jbf.eval_minus(qp.p);

                    mass(gi,gj) += qp.w * jbfval.dot(ibfval);
                }

                edvec3 f{1.0 - qp.p.y(), 1.0 - qp.p.x(), 0};
                rhs(gi) += qp.w * f.dot(ibfval);
            }
        }
    }

    ddvector sol = mass.lu().solve(rhs);

    H5Easy::File file("proj.h5", H5Easy::File::Truncate);
    file.createDataSet("proj/Z", mass);
    file.createDataSet("proj/b", rhs);
    file.createDataSet("proj/x", sol);

    return sol;
}

int
num_shared_vertices(const mesh& msh, const size_t eia, const size_t eib)
{
    assert(eia < msh.edges.size());
    assert(eib < msh.edges.size());
    auto ea = msh.edges[eia];
    auto eb = msh.edges[eib];

    if ( (ea.iv0 == eb.iv0) and (ea.iv1 == eb.iv1) ) {
        return 2;
    }

    if ( (ea.iv0 == eb.iv0) or (ea.iv0 == eb.iv1) or
         (ea.iv1 == eb.iv0) or (ea.iv1 == eb.iv1) ) {
        return 1;
    }

    return 0;
}

void
compute_matrix(const mesh& msh, const std::vector<basis_function>& bfs,
    zdmatrix& Z, const config& cfg)
{
    double freq = cfg.frequency;
    double omega = 2.0*M_PI*freq;
    double k = omega*std::sqrt(MU0*EPS0);
    double inv_ksq = 1./(k*k);

    for (const auto& ibf : bfs) {
        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];
        auto iqps_minus = integrate(msh, iTminus, cfg.degree);
        auto iqps_plus = integrate(msh, iTplus, cfg.degree);

        for (const auto& jbf : bfs) {
            const auto& jTminus = msh.triangles[jbf.itminus];
            const auto& jTplus = msh.triangles[jbf.itplus];
            bool split = num_shared_vertices(msh, ibf.edge_index, jbf.edge_index) != 0;

            auto jqps_minus = split?
                integrate_subtri(msh, jTminus, cfg.degree) :
                integrate(msh, jTminus, cfg.degree);

            auto jqps_plus = split?
                integrate_subtri(msh, jTplus, cfg.degree) :
                integrate(msh, jTplus, cfg.degree);

            std::complex<double> entry = 0.0;
            double prodl = 0.25 * ibf.length * jbf.length * M_1_PI;

            /* iT-, jT- */
            auto inv_AmAm = 1. / (ibf.Aminus * jbf.Aminus);
            for (const auto& iqp : iqps_minus) {
                for (const auto& jqp : jqps_minus) {
                    double v = 0.25*dot( ibf.rho_minus(iqp.p), jbf.rho_minus(jqp.p) );
                    double s = inv_ksq;
                    double w = iqp.w * jqp.w * inv_AmAm;
                    double Rij = norm(iqp.p - jqp.p);
                    std::complex<double> exponent{0, -k*Rij};
                    std::complex<double> g = std::exp(exponent)/Rij;
                    entry += prodl * w * (v - s) * g;
                }
            }

            /* iT-, jT+ */
            auto inv_AmAp = 1. / (ibf.Aminus * jbf.Aplus);
            for (const auto& iqp : iqps_minus) {
                for (const auto& jqp : jqps_plus) {
                    double v = 0.25*dot( ibf.rho_minus(iqp.p), jbf.rho_plus(jqp.p) );
                    double s = inv_ksq;
                    double w = iqp.w * jqp.w * inv_AmAp;
                    double Rij = norm(iqp.p - jqp.p);
                    std::complex<double> exponent{0, -k*Rij};
                    std::complex<double> g = std::exp(exponent)/Rij;
                    entry -= prodl * w * (v - s) * g;
                }
            }

            /* iT+, jT- */
            auto inv_ApAm = 1. / (ibf.Aplus * jbf.Aminus);
            for (const auto& iqp : iqps_plus) {
                for (const auto& jqp : jqps_minus) {
                    double v = 0.25*dot( ibf.rho_plus(iqp.p), jbf.rho_minus(jqp.p) );
                    double s = inv_ksq;
                    double w = iqp.w * jqp.w * inv_ApAm;
                    double Rij = norm(iqp.p - jqp.p);
                    std::complex<double> exponent{0, -k*Rij};
                    std::complex<double> g = std::exp(exponent)/Rij;
                    entry -= prodl * w * (v - s) * g;
                }
            }

            /* iT+, jT+ */
            auto inv_ApAp = 1. / (ibf.Aplus * jbf.Aplus);
            for (const auto& iqp : iqps_plus) {
                for (const auto& jqp : jqps_plus) {
                    double v = 0.25*dot( ibf.rho_plus(iqp.p), jbf.rho_plus(jqp.p) );
                    double s = inv_ksq;
                    double w = iqp.w * jqp.w * inv_ApAp;
                    double Rij = norm(iqp.p - jqp.p);
                    std::complex<double> exponent{0, -k*Rij};
                    std::complex<double> g = std::exp(exponent)/Rij;
                    entry += prodl * w * (v - s) * g;
                }
            }

            Z(jbf.matrix_index, ibf.matrix_index) = entry;
        } // for jbf
    } // for ibf
}

void
compute_rhs(const mesh& msh, const std::vector<basis_function>& bfs,
    zdvector& b, const config& cfg)
{
    double freq = cfg.frequency;
    double omega = 2.0*M_PI*freq;
    double wn = omega*std::sqrt(MU0*EPS0);

    for (const auto& ibf : bfs) {
    
        if (not ibf.interface) {
            continue;
        }
    
        auto itf_tag = ibf.interface.value();
        std::cout << "rhs: itf tag " << itf_tag << " idx " << ibf.matrix_index << "\n";

        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];

        std::complex<double> entry{0.0, -ibf.length/(omega*MU0)};

        b(ibf.matrix_index) = entry;
    }
}

} //namespace mommy

int main(int argc, char **argv)
{
    //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    const char *gmsh_geo_path = nullptr;
    const char *silo_path = "test.silo";
    mommy::config cfg;
    cfg.frequency = 0;
    cfg.degree = 2;

    int opt;
    while ((opt = getopt(argc, argv, "f:g:k:s:")) != -1) {
        switch (opt) {
        case 'f':
            cfg.frequency = std::stod(optarg);
            break;
        case 'g':
            gmsh_geo_path = optarg;
            break;
        case 'k':
            cfg.degree = std::stoull(optarg);
            break;
        case 's':
            silo_path = optarg;
            break;
        default:
            std::cerr << "Invalid argument\n";
            return EXIT_FAILURE;
        }
    }

    if (gmsh_geo_path == nullptr) {
        std::cerr << "No geometry specified (-g)\n";
        return EXIT_FAILURE;
    }

    try {
        gmsh::initialize();
        gmsh::open(gmsh_geo_path);
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

    mommy::silo db(silo_path);
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

    mommy::ddvector vals = mommy::ddvector::Zero(msh.vertices.size());
    for (size_t i = 0; i < msh.vertices.size(); i++) {
        const auto& vtx = msh.vertices[i];
        auto val = std::sin(M_PI*vtx.x())*std::sin(M_PI*vtx.y());
        vals(i) = val;
    }
    db.add_variable("mesh", "vals", vals, mommy::var_centering::nodal);

    mommy::ddfield norms = mommy::ddfield::Zero(msh.triangles.size(), 3);
    for (size_t i = 0; i < msh.triangles.size(); i++) {
        norms.row(i) = normal(msh, msh.triangles[i]);
    }
    db.add_variable("mesh", "normals", norms, mommy::var_centering::zonal);

    std::cout << "Vertices: " << msh.vertices.size() << std::endl;
    std::cout << "Edges:    " << msh.edges.size() << std::endl;
    std::cout << "Cells:    " << msh.triangles.size() << std::endl;
    std::cout << "IntEdges: " << mommy::num_internal_edges(msh) << std::endl;

    std::vector<mommy::basis_function> bfs;
    mommy::make_function_space(msh, bfs);

    #if 0
    for (size_t itri = 0; itri < msh.tbis.size(); itri++) {
        const auto& tbi = msh.tbis[itri];
        std::cout << "Triangle " << itri << "\n";
        for (size_t i = 0; i < 3; i++) {
            std::cout << "Edge " << i << ": ";
            if (tbi[i]) {
                auto bi = *tbi[i];
                std::cout << bi.edge_index << ", ";
                std::cout << bi.sign << ", " << bi.p; 
            }
            std::cout << "\n";
        }
    }
    #endif

    #if 0
    for (size_t i = 0; i < bfs.size(); i++)
    {
        const auto& bf = bfs[i];

        auto Tminus = msh.triangles[bf.itminus];
        auto Tplus  = msh.triangles[bf.itplus];

        std::string fname = "debug/basis_" + std::to_string(bf.edge_index) + "_minus.dat";
        std::ofstream ofs_minus(fname);
        auto qpsminus = mommy::integrate(msh, Tminus, 8);
        for (const auto& qp : qpsminus) {
            auto rho = bf.eval_minus(qp.p);
            ofs_minus << qp.p.x() << " " << qp.p.y() << " " << rho(0) << " " << rho(1) << "\n";
        }

        fname = "debug/basis_" + std::to_string(bf.edge_index) + "_plus.dat";
        std::ofstream ofs_plus(fname);
        auto qpsplus = mommy::integrate(msh, Tplus, 8);
        for (const auto& qp : qpsplus) {
            auto rho = bf.eval_plus(qp.p);
            ofs_plus << qp.p.x() << " " << qp.p.y() << " " << rho(0) << " " << rho(1) << "\n";
        }
    }
    #endif

    auto system_size = mommy::num_internal_edges(msh);
    
    mommy::zdmatrix Z = mommy::zdmatrix::Zero(system_size, system_size);
    mommy::zdvector b = mommy::zdvector::Zero(system_size);

    //#if 0
    std::cout << "Assemblying linear system...\n";
    const auto asm_start{std::chrono::steady_clock::now()};
    mommy::compute_matrix(msh, bfs, Z, cfg);
    mommy::compute_rhs(msh, bfs, b, cfg);
    const auto asm_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
    std::cout << "ASM time: " << asm_elapsed_seconds << " seconds\n";

    std::cout << "Solving linear system...\n";
    const auto start{std::chrono::steady_clock::now()};
    mommy::zdvector x = Z.lu().solve(b);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << "Solve time: " << elapsed_seconds << " seconds\n";

    H5Easy::File file("mommy.h5", H5Easy::File::Truncate);
    file.createDataSet("mommy/Z", Z);
    file.createDataSet("mommy/b", b);
    file.createDataSet("mommy/x", x);
    //#endif

    //mommy::ddvector x = prjtest(msh, bfs);

    std::vector<double> data( msh.triangles.size() );

    mommy::ddfield vdata = mommy::ddfield::Zero( msh.triangles.size(), 3 );

    Eigen::Matrix<double, 3, 1> temp;

    for (const auto& ibf : bfs) {
    
        const auto& Tminus = msh.triangles[ibf.itminus];
        const auto& Tplus = msh.triangles[ibf.itplus];

        auto bar_Tminus = barycenter(msh, Tminus);
        auto bar_Tplus = barycenter(msh, Tplus);

        auto cminus = (ibf.eval_minus(bar_Tminus)*x(ibf.matrix_index)).eval();
        for (size_t i = 0; i < 3; i++) {
            temp(i) = std::real(cminus(i));
        }
        vdata.row(ibf.itminus) += temp;

        auto cplus = (ibf.eval_plus(bar_Tplus)*x(ibf.matrix_index)).eval();
        for (size_t i = 0; i < 3; i++) {
            temp(i) = std::real(cplus(i));
        }
        vdata.row(ibf.itplus) += temp;

        //data[ibf.itminus] += cminus;//cminus.dot(cminus).real();
        //data[ibf.itplus] += cplus;//cplus.dot(cplus).real();
    }

    db.add_variable("mesh", "mag", data, mommy::var_centering::zonal);
    db.add_variable("mesh", "J", vdata, mommy::var_centering::zonal);

    return 0;
}

