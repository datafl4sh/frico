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

#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include "eigen.h"

#include <highfive/H5Easy.hpp>

#include "gmsh.h"

#include "geom_mesh.h"
#include "input_gmsh.h"
#include "output_silo.h"
#include "quadratures.h"
#include "rwg_basis.h"
#include "utils.h"

#define MU0     1.256637061435917e-06
#define EPS0    8.8541878188e-12
namespace frico {

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
                edvec3 ibfval = (bi.sign > 0)?
                    ibf.eval_plus(qp.p) : ibf.eval_minus(qp.p);

                for (size_t j = 0; j < 3; j++) {
                    if (not tbi[j]) { continue; }
                    auto bj = *tbi[j];
                    auto gj = cmap[ bj.edge_index ];
                    const auto& jbf = bfs[ gj ];
                    edvec3 jbfval = (bj.sign > 0)?
                        jbf.eval_plus(qp.p) : jbf.eval_minus(qp.p);

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

void
compute_matrix(const mesh& msh, const std::vector<basis_function>& bfs,
    zdmatrix& Z, const config& cfg)
{
    double freq = cfg.frequency;
    double omega = 2.0*M_PI*freq;
    double k = omega*std::sqrt(MU0*EPS0);
    double inv_ksq = 1./(omega*omega*MU0*EPS0);

    for (const auto& ibf : bfs) {
        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];
        auto iqps_minus = integrate(msh, iTminus, cfg.degree);
        auto iqps_plus = integrate(msh, iTplus, cfg.degree);

        for (const auto& jbf : bfs) {
            const auto& jTminus = msh.triangles[jbf.itminus];
            const auto& jTplus = msh.triangles[jbf.itplus];

            bool split_minus =
                (jbf.itminus == ibf.itminus) or (jbf.itminus == ibf.itplus);
            bool split_plus =
                (jbf.itplus == ibf.itminus) or (jbf.itplus == ibf.itplus);

            auto jqps_minus = split_minus?
                integrate_subtri(msh, jTminus, cfg.degree) :
                integrate(msh, jTminus, cfg.degree);

            auto jqps_plus = split_plus?
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
compute_matrix_approx(const mesh& msh, const std::vector<basis_function>& bfs,
    zdmatrix& Z, const config& cfg)
{
    double freq = cfg.frequency;
    double omega = 2.0*M_PI*freq;
    double k = omega*std::sqrt(MU0*EPS0);
    double inv_ksq = 1./(omega*omega*MU0*EPS0);

    double sqrt_4pi = std::sqrt(4*M_PI);
 
    double dmax = 0.0;

    for (const auto& ibf : bfs) {
        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];
        auto iTminusbar = barycenter(msh, iTminus);
        auto iTplusbar = barycenter(msh, iTplus);

        for (const auto& jbf : bfs) {
            const auto& jTminus = msh.triangles[jbf.itminus];
            const auto& jTplus = msh.triangles[jbf.itplus];
            auto jTminusbar = barycenter(msh, jTminus);
            auto jTplusbar = barycenter(msh, jTplus);

            bool selfmm = (ibf.itminus == jbf.itminus);
            bool selfmp = (ibf.itminus == jbf.itplus);
            bool selfpm = (ibf.itplus == jbf.itminus);
            bool selfpp = (ibf.itplus == jbf.itplus);

            auto sqrt_Am = std::sqrt(measure(msh, iTminus));
            auto sqrt_Ap = std::sqrt(measure(msh, iTplus));
            std::complex jk{0., k};

            std::complex<double> entry = 0.0;

            double rho_cmm = dot( ibf.rho_minus(iTminusbar), jbf.rho_minus(jTminusbar) );
            double Rmm = norm(iTminusbar - jTminusbar);
            std::complex<double> expmm{0, -k*Rmm};
            std::complex<double> Gmm = std::exp(expmm)/(4.0*M_PI*Rmm);
            if (selfmm) {
                Gmm = 0.25*(sqrt_4pi/sqrt_Am - jk) * M_1_PI;
            }
            entry += ibf.length*jbf.length*(0.25*rho_cmm*Gmm - inv_ksq*Gmm);

            double rho_cmp = dot( ibf.rho_minus(iTminusbar), jbf.rho_plus(jTplusbar) );
            double Rmp = norm(iTminusbar - jTplusbar);
            std::complex<double> expmp{0, -k*Rmp};
            std::complex<double> Gmp = std::exp(expmp)/(4.0*M_PI*Rmp);
            if (selfmp) {
                Gmp = 0.25*(sqrt_4pi/sqrt_Am - jk) * M_1_PI;
            }
            entry -= ibf.length*jbf.length*(0.25*rho_cmp*Gmp - inv_ksq*Gmp);

            double rho_cpm = dot( ibf.rho_plus(iTplusbar), jbf.rho_minus(jTminusbar) );
            double Rpm = norm(iTplusbar - jTminusbar);
            std::complex<double> exppm{0, -k*Rpm};
            std::complex<double> Gpm = std::exp(exppm)/(4.0*M_PI*Rpm);
            if (selfpm) {
                Gpm = 0.25*(sqrt_4pi/sqrt_Ap - jk) * M_1_PI;
            }
            entry -= ibf.length*jbf.length*(0.25*rho_cpm*Gpm - inv_ksq*Gpm);

            double rho_cpp = dot( ibf.rho_plus(iTplusbar), jbf.rho_plus(jTplusbar) );
            double Rpp = norm(iTplusbar - jTplusbar);
            std::complex<double> exppp{0, -k*Rpp};
            std::complex<double> Gpp = std::exp(exppp)/(4.0*M_PI*Rpp);
            if (selfpp) {
                Gpp = 0.25*(sqrt_4pi/sqrt_Ap - jk) * M_1_PI;
            }
            entry += ibf.length*jbf.length*(0.25*rho_cpp*Gpp - inv_ksq*Gpp);

            Z(jbf.matrix_index, ibf.matrix_index) = entry;
        } // for jbf
    } // for ibf

    std::cout<< "dmax = "<< dmax << "\n";
}

void
compute_rhs(const mesh& msh, const std::vector<basis_function>& bfs,
    zdvector& b, const config& cfg)
{
    double freq = cfg.frequency;
    double omega = 2.0*M_PI*freq;
    double wn = omega*std::sqrt(MU0*EPS0);

    edvec3 E = {1.0, 0.0, 0.0};

    for (const auto& ibf : bfs) {
    
        if (not ibf.interface) {
            continue;
        }
    
        auto itf_tag = ibf.interface.value();
        std::cout << "rhs: itf tag " << itf_tag << " idx " << ibf.matrix_index << "\n";

        double val = 0.0;

        const auto& Tminus = msh.triangles[ibf.itminus];
        const auto& qpsminus = frico::integrate(msh, Tminus, 1);
        for (auto& qp : qpsminus) {
            val += qp.w * E.dot(ibf.eval_minus(qp.p)); 
        }


        const auto& Tplus = msh.triangles[ibf.itplus];
        const auto& qpsplus = frico::integrate(msh, Tplus, 1);
        for (auto& qp : qpsplus) {
            val += qp.w * E.dot(ibf.eval_plus(qp.p)); 
        }

        auto e = msh.edges[ibf.edge_index];
        auto bar = barycenter(msh, e);

        std::cout << val << " " << ibf.rho_plus(bar) << std::endl;

        std::complex<double> entry{0.0, -ibf.length/(omega*MU0)};

        b(ibf.matrix_index) = entry;//std::complex<double>{0.0, -val/(omega*MU0)};
    }
}

} //namespace frico



bool make_sampling_sphere(frico::mesh& msh, const frico::point& center, double r, double h)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add("sampling");

    gmsh::model::occ::addSphere(center.x(), center.y(), center.z(), r);

    gmsh::model::occ::synchronize();

    gmsh::vectorpair vp;
    gmsh::model::getEntities(vp);
    gmsh::model::mesh::setSize(vp, h);
    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::setOrder(1);

    frico::load_mesh_from_gmsh(msh);

    gmsh::clear();
    gmsh::finalize();
    return true;
}

bool make_sampling_rectangle(frico::mesh& msh, const frico::point& center, double h)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add("sampling");

    //gmsh::model::occ::addCircle(center.x(), center.y(), center.z(), r);
    gmsh::model::occ::addRectangle(-5, -5, 0.0, 10, 10);

    gmsh::model::occ::synchronize();

    gmsh::vectorpair vp;
    gmsh::model::getEntities(vp);
    gmsh::model::mesh::setSize(vp, h);

    gmsh::model::mesh::generate(2);
    //gmsh::model::mesh::setOrder(1);

    frico::load_mesh_from_gmsh(msh);

    gmsh::clear();
    gmsh::finalize();
    return true;
}

int main(int argc, char **argv)
{
    //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    const char *gmsh_geo_path = nullptr;
    const char *silo_path = "test.silo";
    const char *range_expr = nullptr;
    frico::config cfg;
    cfg.frequency = -1;
    cfg.degree = 2;
    bool force_symmetry = false;
    bool approx_matrix = false;

    int opt;
    while ((opt = getopt(argc, argv, "Af:g:k:s:SR:")) != -1) {
        switch (opt) {
        case 'A':
            approx_matrix = true;
            break;
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
        case 'S':
            force_symmetry = true;
            break;
        case 'R':
            range_expr = optarg;
            break;

        default:
            std::cerr << "Invalid argument\n";
            return EXIT_FAILURE;
        }
    }

    /* (1) Check if geometry was specified */
    if (gmsh_geo_path == nullptr) {
        std::cerr << "No geometry specified (-g)\n";
        return EXIT_FAILURE;
    }

    /* (2) Check if frequency was specified, either single or sweep */
    if ( (cfg.frequency <= 0) and (range_expr == nullptr) ) {
        std::cerr << "Simulation frequency not specified (-f or -R)" << std::endl;
        return EXIT_FAILURE;
    }

    /* Generate & load mesh */
    frico::mesh msh;
    if ( not frico::load_mesh_from_gmsh(gmsh_geo_path, msh) ) {
        std::cerr << "Can't load geometry, exiting" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Mesh information: " << std::endl;
    std::cout << "        Vertices: " << msh.vertices.size() << std::endl;
    std::cout << "           Edges: " << msh.edges.size() << std::endl;
    std::cout << "           Cells: " << msh.triangles.size() << std::endl;
    std::cout << "  Internal edges: " << frico::num_internal_edges(msh) << "\n";

    frico::silo db(silo_path);
    db.add_mesh("mesh", msh);

    std::vector<frico::basis_function> bfs;
    frico::make_function_space(msh, bfs);

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
        auto qpsminus = frico::integrate(msh, Tminus, 8);
        for (const auto& qp : qpsminus) {
            auto rho = bf.eval_minus(qp.p);
            ofs_minus << qp.p.x() << " " << qp.p.y() << " " << rho(0) << " " << rho(1) << "\n";
        }

        fname = "debug/basis_" + std::to_string(bf.edge_index) + "_plus.dat";
        std::ofstream ofs_plus(fname);
        auto qpsplus = frico::integrate(msh, Tplus, 8);
        for (const auto& qp : qpsplus) {
            auto rho = bf.eval_plus(qp.p);
            ofs_plus << qp.p.x() << " " << qp.p.y() << " " << rho(0) << " " << rho(1) << "\n";
        }
    }
    #endif


    if (range_expr) {
        auto exp_range = frico::parse_frequency_range(range_expr);
        if ( not exp_range.has_value() ) {
            switch ( exp_range.error() ) {
            case frico::parse_error::invalid_input:
                std::cerr << "Malformed frequency range specification\n";
                return EXIT_FAILURE;
                break;
            case frico::parse_error::out_of_range:
                std::cerr << "Invalid range specification\n";
                return EXIT_FAILURE;
                break;
            default:
                std::unreachable();
            }
        }

        auto range = *exp_range;

        std::ofstream ofs("zplot.txt");

        for (double freq = range.start; freq <= range.end; freq += range.step) {
            cfg.frequency = freq;
            auto system_size = frico::num_internal_edges(msh);
            frico::zdmatrix Z = frico::zdmatrix::Zero(system_size, system_size);
            frico::zdvector b = frico::zdvector::Zero(system_size);

            std::cout << "Assemblying linear system...\n";
            const auto asm_start{std::chrono::steady_clock::now()};
            if (approx_matrix) {
                frico::compute_matrix_approx(msh, bfs, Z, cfg);
            } else {
                frico::compute_matrix(msh, bfs, Z, cfg);
            }
            frico::compute_rhs(msh, bfs, b, cfg);
            const auto asm_end{std::chrono::steady_clock::now()};
            const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
            std::cout << "ASM time: " << asm_elapsed_seconds << " seconds\n";
            if (force_symmetry) {
                Z = (Z+Z.transpose())/2.0;
            }

            std::cout << "Solving linear system...\n";
            const auto start{std::chrono::steady_clock::now()};
            frico::zdvector x = Z.lu().solve(b);
            const auto end{std::chrono::steady_clock::now()};
            const std::chrono::duration<double> elapsed_seconds{end - start};
            std::cout << "Solve time: " << elapsed_seconds << " seconds\n";

            std::complex<double> I = 0.0;
            for (const auto& ibf : bfs) {
                if (not ibf.interface) {
                    continue;
                }
        
                I += ibf.length*x(ibf.matrix_index);
            }
            auto z = 1./I;

            std::complex<double> gamma = (z - 50.0)/(z + 50.0);
            double swr = (1.0 + std::abs(gamma))/(1.0 - std::abs(gamma));
            
            ofs << freq << " " << z.real() << " " << z.imag() << " " << swr << std::endl;

        }

        return 0;
    }


    auto system_size = frico::num_internal_edges(msh);
    
    frico::zdmatrix Z = frico::zdmatrix::Zero(system_size, system_size);
    frico::zdvector b = frico::zdvector::Zero(system_size);

    std::cout << "Assemblying linear system...\n";
    const auto asm_start{std::chrono::steady_clock::now()};
    if (approx_matrix) {
        frico::compute_matrix_approx(msh, bfs, Z, cfg);
    } else {
        frico::compute_matrix(msh, bfs, Z, cfg);
    }
    frico::compute_rhs(msh, bfs, b, cfg);
    const auto asm_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
    std::cout << "ASM time: " << asm_elapsed_seconds << " seconds\n";
    if (force_symmetry) {
        Z = (Z+Z.transpose())/2.0;
    }
  
    H5Easy::File file("frico.h5", H5Easy::File::Truncate);
    file.createDataSet("frico/Z", Z);
    file.createDataSet("frico/b", b);

    std::cout << "Solving linear system...\n";
    const auto start{std::chrono::steady_clock::now()};
    frico::zdvector x = Z.lu().solve(b);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << "Solve time: " << elapsed_seconds << " seconds\n";


    file.createDataSet("frico/x", x);


    std::vector<double> data( msh.triangles.size() );
    std::vector<double> datab( msh.triangles.size() );

    frico::zdfield tri_AJ = frico::zdfield::Zero(msh.triangles.size(), 3);
    frico::zdvector tri_AdivJ = frico::zdvector::Zero(msh.triangles.size());
    frico::zdfield tri_J = frico::ddfield::Zero( msh.triangles.size(), 3 );

    Eigen::Matrix<double, 3, 1> temp;

    for (const auto& ibf : bfs) {
    
        const auto& Tminus = msh.triangles[ibf.itminus];
        const auto& Tplus = msh.triangles[ibf.itplus];

        auto bar_Tminus = barycenter(msh, Tminus);
        auto bar_Tplus = barycenter(msh, Tplus);

        auto A_Tminus = frico::measure(msh, Tminus);
        auto A_Tplus = frico::measure(msh, Tplus);

        frico::ezvec3 Jminus = ibf.eval_minus(bar_Tminus)*x(ibf.matrix_index);
        frico::ezvec3 Jplus = ibf.eval_plus(bar_Tplus)*x(ibf.matrix_index);

        tri_AJ.row(ibf.itminus) += A_Tminus * Jminus;
        tri_AJ.row(ibf.itplus) += A_Tplus * Jplus;

        tri_AdivJ(ibf.itminus) += A_Tminus * ibf.div_minus(bar_Tminus)*x(ibf.matrix_index);
        tri_AdivJ(ibf.itplus) += A_Tplus * ibf.div_plus(bar_Tplus)*x(ibf.matrix_index);

        tri_J.row(ibf.itminus) += Jminus;
        tri_J.row(ibf.itplus) += Jplus;
    }

    db.add_variable("mesh", "J", tri_J, frico::var_centering::zonal);

    frico::zdfield Exy = frico::zdfield::Zero(360, 3);
    frico::zdfield Eyz = frico::zdfield::Zero(360, 3);
    frico::zdfield Exz = frico::zdfield::Zero(360, 3);
    

    double freq = cfg.frequency;
    double omega = 2.0*M_PI*freq;
    double k = omega*std::sqrt(MU0*EPS0);
    std::complex<double> jomega{0.0, omega};



    for (int deg = 0; deg < 360; deg++) {
        double theta = deg*M_PI/180;
        double R = 5.0;
        
        for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
            const auto& tri = msh.triangles[itri];
            auto bar = barycenter(msh, tri);
            frico::ezvec3 J = tri_AJ.row(itri);
            std::complex<double> divJ = tri_AdivJ(itri);

            /* XY */ {
                frico::point Pxy{ R*std::cos(theta), R*std::sin(theta), 0.0 };
                frico::vec3 vR = Pxy - bar;
                double R = norm(vR);
                std::complex<double> jkR{0.0, k*R};
                std::complex<double> g = (std::exp(-jkR)/(4*M_PI*R));
                Exy.row(deg) += -jomega*MU0*(J - divJ*(1.0+jkR)*vR.to_eigen()/(k*k*R*R))*g;
            }

            /* YZ */ {
                frico::point Pyz{ 0.0, R*std::cos(theta), R*std::sin(theta) };
                frico::vec3 vR = Pyz - bar;
                double R = norm(vR);
                std::complex<double> jkR{0.0, k*R};
                std::complex<double> g = (std::exp(-jkR)/(4*M_PI*R));
                Eyz.row(deg) += -jomega*MU0*(J - divJ*(1.0+jkR)*vR.to_eigen()/(k*k*R*R))*g;
            }

            /* XZ */ {
                frico::point Pxz{ R*std::cos(theta), 0.0, -R*std::sin(theta) };
                frico::vec3 vR = Pxz - bar;
                double R = norm(vR);
                std::complex<double> jkR{0.0, k*R};
                std::complex<double> g = (std::exp(-jkR)/(4*M_PI*R));
                Exz.row(deg) += -jomega*MU0*(J - divJ*(1.0+jkR)*vR.to_eigen()/(k*k*R*R))*g;
            }
        }
    }

    std::ofstream ofs_nec("neccomp.csv");
    ofs_nec << "X;Y;Z;EXM;EXP;EYM;EYP;EZM;EZP" << std::endl;
    double sx = -5.0;
    double sz = -5.0;
    double dx = 0.25;
    double dz = 0.25;
    for (int i = 0; i < 41; i++) {
        for (int j = 0; j < 41; j++) {
            double x = sx + j*dx;
            double z = sz + i*dz;
            frico::point pt{x, 0.0, z};

            frico::ezvec3 locE = frico::ezvec3::Zero();
            for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
                const auto& tri = msh.triangles[itri];
                auto bar = barycenter(msh, tri);
                frico::ezvec3 J = tri_AJ.row(itri);
                std::complex<double> divJ = tri_AdivJ(itri);

                frico::vec3 vR = pt - bar;
                double R = norm(vR);
                std::complex<double> jkR{0.0, k*R};
                std::complex<double> g = (std::exp(-jkR)/(4*M_PI*R));
                locE += -jomega*MU0*(J - divJ*(1.0+jkR)*vR.to_eigen()/(k*k*R*R))*g;
            }
            
            ofs_nec << pt.x() << ";" << pt.y() << ";" << pt.z() << ";";
            ofs_nec << std::abs(locE(0)) << ";" << std::arg(locE(0)) << ";";
            ofs_nec << std::abs(locE(1)) << ";" << std::arg(locE(1)) << ";";
            ofs_nec << std::abs(locE(2)) << ";" << std::arg(locE(2)) << "\n";
        }
    }
    
    double magExy_max = 0.0;
    double magEyz_max = 0.0;
    double magExz_max = 0.0;
    for (int i = 0; i < 360; i++) {
        { double magExy = std::sqrt(std::real(Exy.row(i).dot(Exy.row(i))));
          magExy_max = std::max(magExy_max, magExy); }
        { double magEyz = std::sqrt(std::real(Eyz.row(i).dot(Eyz.row(i))));
          magEyz_max = std::max(magEyz_max, magEyz); }
        { double magExz = std::sqrt(std::real(Exz.row(i).dot(Exz.row(i))));
          magExz_max = std::max(magExz_max, magExz); }
    }

    std::ofstream ofs("polar.txt");
    for (int i = 0; i < 360; i++) {
        double magExy = std::sqrt(std::real(Exy.row(i).dot(Exy.row(i))));
        double magEyz = std::sqrt(std::real(Eyz.row(i).dot(Eyz.row(i))));
        double magExz = std::sqrt(std::real(Exz.row(i).dot(Exz.row(i))));
        ofs << i << " " << magExy << " " << 20*std::log10(magExy/magExy_max)
                 << " " << magEyz << " " << 20*std::log10(magEyz/magEyz_max)
                 << " " << magExz << " " << 20*std::log10(magExz/magExz_max)
                 << std::endl;
    }

    #if 0
    frico::mesh samplingsphere;
    //make_sampling_sphere(samplingsphere, {0,0,0}, 0.75, 0.05);
    make_sampling_rectangle(samplingsphere, {0,0,0}, 0.025);
    auto nspoints = samplingsphere.vertices.size();
    std::vector<double> sEmag( nspoints );
    frico::zdfield vEmag = frico::zdfield::Zero( nspoints, 3 );
    std::cout << "postpro begin\n";
    for (size_t i = 0; i < nspoints; i++) {
        const auto& spt = samplingsphere.vertices[i]; 
        frico::ezvec3 locE = frico::ezvec3::Zero();
        for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
            const auto& tri = msh.triangles[itri];
            frico::vec3 bar = frico::barycenter(msh, tri);
            double R = norm(spt - bar);
            std::complex<double> jkR{0, k*R};
            std::complex<double> g = (std::exp(-jkR)/(4*M_PI*R));
            frico::ezvec3 J = tri_AJ.row(itri);
            std::complex<double> divJ = tri_AdivJ(itri);
            locE += -jomega*MU0*(J - divJ*(1.0+jkR)*(spt-bar).to_eigen()/(k*k*R*R))*g;
        }
        sEmag[i] = std::sqrt(std::real(locE.dot(locE)));
        vEmag.row(i) = locE;
    }
    std::cout << "postpro end\n";

    db.add_mesh("sampling", samplingsphere);
    db.add_variable("sampling", "magE", sEmag, frico::var_centering::nodal);
    db.add_variable("sampling", "E", vEmag, frico::var_centering::nodal);
    
    
    for (int i = 0; i < 100; i++) {
        std::string stepfname = "step_";
        stepfname = stepfname + std::to_string(i) + ".silo";
        frico::silo db2(stepfname);

        auto dT = i*(1./cfg.frequency)/100.0;

        frico::ddfield E = (vEmag*std::exp(jomega*dT)).real();

        db2.add_mesh("sampling", samplingsphere);
        db2.add_variable("sampling", "E", E, frico::var_centering::nodal);
    }
    #endif

    std::complex<double> I = 0.0;
    for (const auto& ibf : bfs) {
    
        if (not ibf.interface) {
            continue;
        }
    
        I += ibf.length*x(ibf.matrix_index);
    }

    auto z = 1./I;
    std::complex<double> gamma = (z - 50.0)/(z + 50.0);
    double swr = (1.0 + std::abs(gamma))/(1.0 - std::abs(gamma));
    std::cout << "Impedance: " << z << std::endl;
    std::cout << "      SWR: " << swr << std::endl;

    return 0;
}

