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

#include <print>
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
#include "sources.h"
#include "constants.h"
#include "bem_maxwell.h"












int main(int argc, char **argv)
{
    //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    #ifndef _OPENMP
    std::print(
"FRICO v0.0 3D MoM solver - Matteo `IV3IWE` Cicuttin (C) 2025-2026\n\n";
    );
    #else
    std::print(
"FRICO v0.0 3D MoM solver - Matteo `IV3IWE` Cicuttin (C) 2025-2026 [OpenMP]\n\n"
    );
    #endif

    frico::maxwell::simulation sim;

    const char *arg_geo_path = nullptr;
    const char *arg_skiptags = nullptr;
    const char *arg_frequency = nullptr;
    const char *arg_range_expr = nullptr;
    const char *arg_simname = "default";

    int opt;
    while ((opt = getopt(argc, argv, "Af:g:k:n:s:SR:x:")) != -1) {
        switch (opt) {
        case 'A':
            sim.cfg.approx_matrix = true;
            break;
        case 'f':
            arg_frequency = optarg;
            break;
        case 'g':
            arg_geo_path = optarg;
            break;
        case 'k':
            sim.cfg.degree = std::stoull(optarg);
            break;
        case 'n':
            arg_simname = optarg;
            break;
        case 's':
            //sim.cfg.silo_path = optarg;
            break;
        case 'S':
            sim.cfg.force_symmetry = true;
            break;
        case 'R':
            arg_range_expr = optarg;
            break;
        case 'x':
            arg_skiptags = optarg;
            break;

        default:
            std::println(stderr, "Invalid argument");
            return EXIT_FAILURE;
        }
    }

    /* (1) Check if geometry was specified */
    if (arg_geo_path == nullptr) {
        std::println(stderr, "No geometry specified (-g)");
        return EXIT_FAILURE;
    }

    /* (2) Check if frequency was specified, either single or sweep */
    auto opt_freqs = frico::parse_frequency_parameters(arg_frequency, arg_range_expr);
    if (not opt_freqs) {
        return EXIT_FAILURE;
    }

    /* (3) If there is a list of surfaces to skip, process it */
    if (arg_skiptags) {
        auto exp_skiptags = frico::parse_integer_list(arg_skiptags);
        if (exp_skiptags.has_value()) {
            sim.skiptags = *exp_skiptags;
        } else {
            std::println(stderr, "Error parsing the argument of -x");
            return EXIT_FAILURE;
        }
    }


    if ( not frico::maxwell::init_simulation(sim, arg_simname, arg_geo_path) ) {
        return EXIT_FAILURE;
    }

    frico::maxwell::init_sweep(sim, *opt_freqs);
    frico::maxwell::run(sim);

    return EXIT_SUCCESS;
}




#if 0

int main(int argc, char **argv)
{
    //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    std::cout <<
        "FRICO v0.0 3D MoM solver - Matteo Cicuttin [IV3IWE] (C) 2025-2026\n\n";
    
    frico::config cfg;
    cfg.frequency = -1;
    cfg.degree = 2;

    const char *arg_skiptags = nullptr;
    const char *arg_range_expr = nullptr;

    int opt;
    while ((opt = getopt(argc, argv, "Af:g:k:s:SR:x:")) != -1) {
        switch (opt) {
        case 'A':
            cfg.approx_matrix = true;
            break;
        case 'f':
            cfg.frequency = std::stod(optarg);
            break;
        case 'g':
            cfg.gmsh_geo_path = optarg;
            break;
        case 'k':
            cfg.degree = std::stoull(optarg);
            break;
        case 's':
            cfg.silo_path = optarg;
            break;
        case 'S':
            cfg.force_symmetry = true;
            break;
        case 'R':
            arg_range_expr = optarg;
            break;
        case 'x':
            arg_skiptags = optarg;
            break;

        default:
            std::cerr << "Invalid argument\n";
            return EXIT_FAILURE;
        }
    }

    /* (1) Check if geometry was specified */
    if (cfg.gmsh_geo_path == nullptr) {
        std::cerr << "No geometry specified (-g)\n";
        return EXIT_FAILURE;
    }

    /* (2) Check if frequency was specified, either single or sweep */
    if ( (cfg.frequency <= 0) and (arg_range_expr == nullptr) ) {
        std::cerr << "Simulation frequency not specified (-f or -R)" << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<int> skiptags;
    if (arg_skiptags) {
        auto exp_skiptags = frico::parse_integer_list(arg_skiptags);
        if (not exp_skiptags) {
            std::cerr << "Error parsing argument of -x option" << std::endl;
            return EXIT_FAILURE;
        }
        skiptags = *exp_skiptags;
    }

    /* Generate & load mesh */
    frico::mesh msh;
    auto meshok = frico::load_from_gmsh(cfg.gmsh_geo_path, msh, skiptags);
    if ( not meshok ) {
        auto err = meshok.error();
        if (err.errtype == frico::meshing_error::gmsh_issue) {
            std::cerr << "Can't load geometry, exiting" << std::endl;
            return EXIT_FAILURE;
        }
        if (err.errtype == frico::meshing_error::bad_connectivity) {
            std::cerr <<
                "The mesh connectivity is not valid because an edge shared by\n"
                "more than two triangles was detected. This does not permit\n"
                "to construct the RWG basis.\n"
                "The tags of the involved surfaces are " << err.tminus_tag
                << ", " << err.tplus_tag << " and " << err.offending_tag << ".\n";
            return EXIT_FAILURE;
        }
    }

    std::cout << "Mesh information: " << std::endl;
    std::cout << "        Vertices: " << msh.vertices.size() << std::endl;
    std::cout << "           Edges: " << msh.edges.size() << std::endl;
    std::cout << "           Cells: " << msh.triangles.size() << std::endl;
    std::cout << "  Internal edges: " << frico::num_internal_edges(msh) << "\n";

    frico::silo db(cfg.silo_path);
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


    if (arg_range_expr) {
        auto exp_range = frico::parse_frequency_range(arg_range_expr);
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
                std::cerr << "Unreachable branch taken" << std::endl;
                std::terminate();
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
            if (cfg.approx_matrix) {
                frico::compute_matrix_approx(msh, bfs, Z, cfg);
            } else {
                frico::compute_matrix(msh, bfs, Z, cfg);
            }
            frico::compute_rhs(msh, bfs, b, cfg);
            const auto asm_end{std::chrono::steady_clock::now()};
            const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
            std::cout << "ASM time: " << asm_elapsed_seconds << " seconds\n";
            if (cfg.force_symmetry) {
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
        
                //if (ibf.interface.value() == 3)
                //    I += ibf.length*x(ibf.matrix_index);
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
    if (cfg.approx_matrix) {
        frico::compute_matrix_approx(msh, bfs, Z, cfg);
    } else {
        frico::compute_matrix(msh, bfs, Z, cfg);
    }
    frico::compute_rhs(msh, bfs, b, cfg);
    const auto asm_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
    std::cout << "ASM time: " << asm_elapsed_seconds << " seconds\n";
    if (cfg.force_symmetry) {
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

#endif