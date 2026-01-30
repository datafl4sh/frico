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

#include <fstream>
#include <chrono>

#include "bem_maxwell.h"
#include "input_gmsh.h"
#include "quadratures.h"
#include "output_silo.h"
#include "constants.h"
#include "utils.h"

#include "gmsh.h"

namespace frico::maxwell {

void
compute_matrix(simulation& sim, size_t ctx_number)
{
    auto& context = sim.contexts[ctx_number];
    const auto& msh = sim.msh;

    double freq = context.frequency;
    double omega = 2.0*M_PI*freq;
    double k = omega*std::sqrt(MU0*EPS0);
    double inv_ksq = 1./(omega*omega*MU0*EPS0);

    #pragma omp parallel for
    for (const auto& ibf : sim.bfuncs) {
        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];
        auto iqps_minus = integrate(msh, iTminus, sim.cfg.degree);
        auto iqps_plus = integrate(msh, iTplus, sim.cfg.degree);

        for (const auto& jbf : sim.bfuncs) {
            const auto& jTminus = msh.triangles[jbf.itminus];
            const auto& jTplus = msh.triangles[jbf.itplus];

            bool split_minus =
                (jbf.itminus == ibf.itminus) or (jbf.itminus == ibf.itplus);
            bool split_plus =
                (jbf.itplus == ibf.itminus) or (jbf.itplus == ibf.itplus);

            auto jqps_minus = split_minus?
                integrate_subtri(msh, jTminus, sim.cfg.degree) :
                integrate(msh, jTminus, sim.cfg.degree);

            auto jqps_plus = split_plus?
                integrate_subtri(msh, jTplus, sim.cfg.degree) :
                integrate(msh, jTplus, sim.cfg.degree);

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

            context.Z(jbf.matrix_index, ibf.matrix_index) = entry;
        } // for jbf
    } // for ibf
}

void
compute_matrix_approx(simulation& sim, size_t ctx_number)
{
    auto& context = sim.contexts[ctx_number];
    const auto& msh = sim.msh;

    double freq = context.frequency;
    double omega = 2.0*M_PI*freq;
    double k = omega*std::sqrt(MU0*EPS0);
    double inv_ksq = 1./(omega*omega*MU0*EPS0);

    double sqrt_4pi = std::sqrt(4*M_PI);

    #pragma omp parallel for
    for (const auto& ibf : sim.bfuncs) {
        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];
        auto iTminusbar = barycenter(msh, iTminus);
        auto iTplusbar = barycenter(msh, iTplus);

        for (const auto& jbf : sim.bfuncs) {
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

            context.Z(jbf.matrix_index, ibf.matrix_index) = entry;
        } // for jbf
    } // for ibf
}

void
compute_rhs(simulation& sim, size_t ctx_number)
{
    auto& context = sim.contexts[ctx_number];

    double freq = context.frequency;
    double omega = 2.0*M_PI*freq;

    for (const auto& ibf : sim.bfuncs) {
    
        if (not ibf.interface) {
            continue;
        }

        //if (ibf.interface.value() == 3) {
            std::complex<double> entry{0.0, -ibf.length/(omega*MU0)};
            context.V(ibf.matrix_index) = entry;
        //}
    }
}

bool init_simulation(simulation& sim, const std::string& name,
    const std::string& geo_path)
{
    sim.name = name;

    auto meshok = frico::load_from_gmsh(geo_path, sim.msh, sim.skiptags);
    if ( not meshok ) {
        auto err = meshok.error();
        if (err.errtype == frico::meshing_error::gmsh_issue) {
            std::cerr << "Can't load geometry, exiting" << std::endl;
            return false;
        }
        if (err.errtype == frico::meshing_error::bad_connectivity) {
            std::cerr <<
                "The mesh connectivity is not valid because an edge shared by\n"
                "more than two triangles was detected. This does not permit\n"
                "to construct the RWG basis.\n"
                "The tags of the involved surfaces are " << err.tminus_tag
                << ", " << err.tplus_tag << " and " << err.offending_tag << ".\n";
            return false;
        }
    }

    std::cout << "Mesh information: " << std::endl;
    std::cout << "        Vertices: " << sim.msh.vertices.size() << std::endl;
    std::cout << "           Edges: " << sim.msh.edges.size() << std::endl;
    std::cout << "           Cells: " << sim.msh.triangles.size() << std::endl;
    std::cout << "  Internal edges: " << num_internal_edges(sim.msh) << "\n";

    make_function_space(sim.msh, sim.bfuncs);

    return true;
}

bool init_sweep(simulation& sim, const frequency_range& freqs)
{
    size_t ctx_number = 0;
    for (double freq = freqs.start; freq <= freqs.end; freq += freqs.step) {
        freq_context context;
        context.ctx_number = ctx_number++;
        context.frequency = freq;
        sim.contexts.push_back(context);
    }

    return true;
}

bool run_context(simulation& sim, size_t ctx_number)
{
    freq_context& context = sim.contexts[ctx_number];

    std::cout << "Sweep frequency " << ctx_number << ": ";
    std::cout << context.frequency << " Hz" << std::endl; 

    auto system_size = num_internal_edges(sim.msh);
    context.Z = zdmatrix::Zero(system_size, system_size);
    context.V = zdvector::Zero(system_size);

    std::cout << "  Assemblying linear system..." << std::flush;
    const auto asm_start{std::chrono::steady_clock::now()};
    if (sim.cfg.approx_matrix) {
        compute_matrix_approx(sim, ctx_number);
    } else {
        compute_matrix(sim, ctx_number);
    }
    
    compute_rhs(sim, ctx_number);
    
    const auto asm_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> asm_elapsed_seconds{asm_end - asm_start};
    std::cout << asm_elapsed_seconds << " seconds\n";
    if (sim.cfg.force_symmetry) {
        context.Z = (context.Z+context.Z.transpose())/2.0;
    }

    std::cout << "  Solving linear system..." << std::flush;
    const auto start{std::chrono::steady_clock::now()};
    context.I = context.Z.lu().solve(context.V);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << elapsed_seconds << " seconds\n";

    context.tri_AJ = zdfield::Zero(sim.msh.triangles.size(), 3);
    context.tri_AdivJ = zdvector::Zero(sim.msh.triangles.size());
    context.tri_J = zdfield::Zero(sim.msh.triangles.size(), 3);

    for (const auto& bf : sim.bfuncs) {
        const auto& Tminus = sim.msh.triangles[bf.itminus];
        const auto& Tplus = sim.msh.triangles[bf.itplus];

        auto bar_Tminus = barycenter(sim.msh, Tminus);
        auto bar_Tplus = barycenter(sim.msh, Tplus);

        auto A_Tminus = measure(sim.msh, Tminus);
        auto A_Tplus = measure(sim.msh, Tplus);

        std::complex<double> Iedge = context.I(bf.matrix_index);

        ezvec3 Jminus = bf.eval_minus(bar_Tminus) * Iedge;
        ezvec3 Jplus = bf.eval_plus(bar_Tplus) * Iedge;

        context.tri_AJ.row(bf.itminus) += A_Tminus * Jminus;
        context.tri_AJ.row(bf.itplus) += A_Tplus * Jplus;

        context.tri_AdivJ(bf.itminus) += A_Tminus * bf.div_minus(bar_Tminus) * Iedge;
        context.tri_AdivJ(bf.itplus) += A_Tplus * bf.div_plus(bar_Tplus) * Iedge;

        context.tri_J.row(bf.itminus) += Jminus;
        context.tri_J.row(bf.itplus) += Jplus;
    }

    std::complex<double> totI = 0.0;
    std::complex<double> totP = 0.0;
    for (const auto& ibf : sim.bfuncs) {
        if (not ibf.interface) {
            continue;
        }

        //if (ibf.interface.value() == 3)
            //auto V = context.V(ibf.matrix_index);
            auto I = ibf.length*context.I(ibf.matrix_index);
            totP += 0.5*1.0*conj(I);
            totI += I;
    }
    auto z = 1./totI;

    context.P = totP;

    std::complex<double> gamma = (z - 50.0)/(z + 50.0);
    double swr = (1.0 + std::abs(gamma))/(1.0 - std::abs(gamma));
            
    std::cout << z << " " << swr << " " << totP << std::endl;

    return true;
}

bool make_sampling_sphere(mesh& msh, const point& center,
    double r, double h)
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

    if ( not load_from_gmsh(msh, load_mode::quick) ) {
        return false;
    }

    gmsh::clear();
    gmsh::finalize();
    return true;
}

bool make_sampling_grid(frico::mesh& msh, const frico::point& c,
    double r, double h)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add("sampling");

    int tagxy = gmsh::model::occ::addRectangle(c.x()-r, c.y()-r, c.z(), 2*r, 2*r);
    int tagyz = gmsh::model::occ::addRectangle(c.x()-r, c.y()-r, c.z(), 2*r, 2*r);
    int tagxz = gmsh::model::occ::addRectangle(c.x()-r, c.y()-r, c.z(), 2*r, 2*r);
    gmsh::model::occ::rotate({{2, tagyz}},
        c.x(), c.y(), c.z(), c.x(), c.y()+1, c.z(), M_PI/2 );
    gmsh::model::occ::rotate({{2, tagxz}},
        c.x(), c.y(), c.z(), c.x()+1, c.y(), c.z(), M_PI/2 );

    gmsh::vectorpair out;
    std::vector<gmsh::vectorpair> outmap;
    gmsh::model::occ::fuse({ {2,tagxy}  }, {{2,tagyz}, {2,tagxz}}, out, outmap);

    gmsh::model::occ::synchronize();

    gmsh::vectorpair vp;
    gmsh::model::getEntities(vp);
    gmsh::model::mesh::setSize(vp, h);

    gmsh::model::mesh::generate(2);
    //gmsh::model::mesh::setOrder(1);

    if ( not load_from_gmsh(msh, load_mode::quick) ) {
        return false;
    }

    gmsh::clear();
    gmsh::finalize();
    return true;
}


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

bool
make_radiation_diagrams(const simulation& sim, size_t ctx_number,
    const mesh& smpmsh, ddvector& gain)
{
    const freq_context& context = sim.contexts[ctx_number];

    gain = ddvector::Zero(smpmsh.vertices.size());

    #pragma omp parallel for
    for (size_t i = 0; i < smpmsh.vertices.size(); i++) {
        const auto& p = smpmsh.vertices[i];
        auto [locE, locH] = eval_fields(sim, ctx_number, p);
        auto R = norm(p);
        ezvec3 S = 0.5*locE.cross(locH.conjugate());
        double Prad = std::real(std::sqrt(S.dot(S)));
        double G = 4*M_PI*R*R*Prad/std::real(context.P);
        gain(i) = G;
    }

    return true;
}

bool make_radiation_diagrams(const simulation& sim, size_t ctx_number,
    const point& center, double radius)
{
    const freq_context& context = sim.contexts[ctx_number];

    ddvector Gxy = ddvector::Zero(360);
    ddvector Gyz = ddvector::Zero(360);
    ddvector Gxz = ddvector::Zero(360);

    double maxG = -std::numeric_limits<double>::infinity();

    #pragma omp parallel for
    for (int deg = 0; deg < 360; deg++) {
        double theta = deg*M_PI/180;
        double c = std::cos(theta);
        double s = std::sin(theta);
        double R = 5.0;
  
        /* XY */ {
            point Pxy{ R*c, R*s, 0.0 };
            auto [locE, locH] = eval_fields(sim, ctx_number, Pxy);
            ezvec3 S = 0.5*locE.cross(locH.conjugate());
            double Prad = std::real(std::sqrt(S.dot(S)));
            double G = 4*M_PI*R*R*Prad/std::real(context.P);
            Gxy(deg) = G;
            maxG = std::max(G, maxG);
        }

        /* YZ */ {
            point Pyz{ 0.0, R*c, R*s };
            auto [locE, locH] = eval_fields(sim, ctx_number, Pyz);
            ezvec3 S = 0.5*locE.cross(locH.conjugate());
            double Prad = std::real(std::sqrt(S.dot(S)));
            double G = 4*M_PI*R*R*Prad/std::real(context.P);
            Gyz(deg) = G;
            maxG = std::max(G, maxG);
        }

        /* XZ */ {
            point Pxz{ -R*c, 0.0, R*s };
            auto [locE, locH] = eval_fields(sim, ctx_number, Pxz);
            ezvec3 S = 0.5*locE.cross(locH.conjugate());
            double Prad = std::real(std::sqrt(S.dot(S)));
            double G = 4*M_PI*R*R*Prad/std::real(context.P);
            Gxz(deg) = G;
            maxG = std::max(G, maxG);
        }
    }

    std::cout << "Max gain: " << 10*std::log10(maxG) << " dB\n";

    std::string filename = "polar_" + std::to_string(ctx_number) + ".txt";
    std::ofstream ofs(filename);
    for (int i = 0; i < 360; i++) {
        double magGxy = Gxy(i);
        double magGyz = Gyz(i);
        double magGxz = Gxz(i);
        ofs << i << " " << magGxy << " " << 10*std::log10(magGxy)
                 << " " << magGyz << " " << 10*std::log10(magGyz)
                 << " " << magGxz << " " << 10*std::log10(magGxz)
                 << std::endl;
    }

    return true;
}

bool eval_fields(const simulation& sim, size_t ctx_number,
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

    return true;
}

bool postpro_context(simulation& sim, size_t ctx_number)
{
    freq_context& context = sim.contexts[ctx_number];

    std::string filename =
        sim.name + "_" + std::to_string(ctx_number) + ".silo";

    silo db;
    db.open(filename);
    db.add_mesh("mesh", sim.msh);
    db.add_variable("mesh", "J", context.tri_J, var_centering::zonal);

    mesh smpmsh;
    zdfield E, H;
    make_sampling_grid(smpmsh, {0,0,0}, 5, 0.1);
    db.add_mesh("sampling", smpmsh);
    eval_fields(sim, ctx_number, smpmsh, E, H);

    zdvector Z = zdvector::Zero(smpmsh.vertices.size());
    for (size_t i = 0; i < smpmsh.vertices.size(); i++) {
        ezvec3 lE = E.row(i);
        ezvec3 lH = H.row(i);
        Z(i) = std::sqrt(lE.dot(lE))/std::sqrt(lH.dot(lH));
    }

    db.add_variable("sampling", "E", E, var_centering::nodal);
    db.add_variable("sampling", "H", H, var_centering::nodal);
    db.add_variable("sampling", "Z", Z, var_centering::nodal);

    mesh smpsph;
    make_sampling_sphere(smpsph, {0.0, 0.0, 0.0}, 5, 0.5);
    ddvector gain;
    make_radiation_diagrams(sim, ctx_number, smpsph, gain);
    db.add_mesh("gainsmp", smpsph);
    db.add_variable("gainsmp", "gain", gain, var_centering::nodal);

    db.close();


    make_radiation_diagrams(sim, ctx_number, {0,0,0}, 5.0);


    
    return true;
}

bool run(simulation& sim)
{
    for (size_t ctx_num = 0; ctx_num < sim.contexts.size(); ctx_num++) {
        run_context(sim, ctx_num);
        postpro_context(sim, ctx_num);

        /* dealloc matrix when we're done */
        sim.contexts[ctx_num].Z.resize(0,0);
    }
    return true;
}

}