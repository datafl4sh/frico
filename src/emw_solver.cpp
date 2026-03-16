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

#include <print>
#include <fstream>
#include <chrono>
#include <set>
#include <list>

#include "eigen.h"
#include <highfive/H5Easy.hpp>

#include "emw_solver.h"
#include "input_gmsh.h"
#include "quadratures.h"
#include "output_silo.h"
#include "constants.h"
#include "utils.h"

#include "gmsh.h"

#include "emw_postpro_common.h"
#include "emw_postpro_delta_gap.h"

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

bool
reorient_deltagap_edges(const simulation& sim, const delta_gap& dg,
    std::vector<double>& signs)
{
    /* Get all the edges involved in the delta-gap */
    bool fix_needed = false;
    std::list<std::pair<edge, size_t>> edges;
    for (const auto& ibf : sim.bfuncs) {
        if (not ibf.interface) {
            continue;
        }
        for (const auto& itf : dg.interfaces) {
            if (*ibf.interface == itf) {
                const auto& Tplus = sim.msh.triangles[ibf.itplus];
                const auto& Tminus = sim.msh.triangles[ibf.itminus];
                if (Tplus.tag <= Tminus.tag) {
                    fix_needed = true;
                }
                edges.push_back(
                    {sim.msh.edges[ibf.edge_index], ibf.matrix_index}
                );
            }
        }
    }

    if ( not (fix_needed or sim.cfg.force_reorient_deltagap) ) {
        return true;
    }

    std::println("Reorientation triggered on delta-gap '{}' with tag {} {}",
        dg.name, dg.phys_entity,
        sim.cfg.force_reorient_deltagap ? "[forced]" : ""
    );

    if (edges.size() == 0) {
        std::println("No edges identified for the specified delta-gap");
        return false;
    }

    struct nodepair {
        size_t first;
        size_t second;
        size_t bf_index;
    };
    std::list<nodepair> seglist;

    /* Init with the first edge of the set */
    auto eitor = edges.begin();
    seglist.push_back({eitor->first.iv0, eitor->first.iv1, eitor->second});
    edges.erase(eitor);

    /* Iterate on the remaining edges and connect them in the
     * right place with the appropriate orientation. This algorithm
     * is quite naive, but for now does the job; we do not expect
     * delta-gaps with hudreds of edges */
    while (edges.size() > 0) {
        auto len_before = edges.size();
        for (eitor = edges.begin(); eitor != edges.end(); /**/) {
            assert(seglist.size() > 0);
            auto f = seglist.front();
            auto b = seglist.back();
            auto cur = eitor++;

            if (edges.size() > 0 and cur->first.iv0 == f.first) {
                seglist.push_front({cur->first.iv1, cur->first.iv0, cur->second});
                edges.erase(cur);
            }
            if (edges.size() > 0 and cur->first.iv1 == f.first) {
                seglist.push_front({cur->first.iv0, cur->first.iv1, cur->second});
                edges.erase(cur);
            }
            if (edges.size() > 0 and cur->first.iv0 == b.second) {
                seglist.push_back({cur->first.iv0, cur->first.iv1, cur->second});
                edges.erase(cur);
            }
            if (edges.size() > 0 and cur->first.iv1 == b.second) {
                seglist.push_back({cur->first.iv1, cur->first.iv0, cur->second});
                edges.erase(cur);
            }
        }
        auto len_after = edges.size();

        if (len_after == len_before) {
            std::println("Invalid source: Entities involved in the "
                "delta-gap are not topologically connected");
            return false;
        }
    }
    
    /* OK, we have an ordered chain of edges forming the delta-gap. Compute
     * the vector from the barycenter of the first edge to the barycenter
     * of its associated T+. This defines the positive direction. */
    signs.resize( sim.bfuncs.size() );
    assert(seglist.size() > 0);
    const auto& seg0 = seglist.front();
    const auto& bf0 = sim.bfuncs[ seg0.bf_index ];
    const auto& Tplus0 = sim.msh.triangles[bf0.itplus];
    const auto& e0 = sim.msh.edges[bf0.edge_index];
    auto dir = barycenter(sim.msh, Tplus0) - barycenter(sim.msh, e0);
    double sign = +1;

    for (auto& seg : seglist) {
        /* Now start iterating on the edge chain and verify that 
         * adjacent basis functions have the same orientation
         * (positive dot product). If not, change sign and
         * record that. */
        const auto& bf = sim.bfuncs[ seg.bf_index ];
        const auto& Tplus = sim.msh.triangles[bf.itplus];
        const auto& e = sim.msh.edges[bf.edge_index];
        auto curdir = barycenter(sim.msh, Tplus) - barycenter(sim.msh, e);

        if (dot(dir, curdir) < 0) {
            sign *= -1;
        }

        signs[bf.matrix_index] = sign;
        dir = curdir;
    }

    if (sim.cfg.verbose) {
        std::print("Edges of the delta-gap and orientation: ");
        for (auto& seg : seglist) {
            std::print("[({} {}) {:+}] ", seg.first, seg.second, sign );
        }
        std::println();
    }

    return true;
}

void
update_rhs(simulation& sim, size_t ctx_number, const delta_gap& dg)
{

    auto& context = sim.contexts[ctx_number];

    double freq = context.frequency;
    double omega = 2.0*M_PI*freq;

    for (const auto& ibf : sim.bfuncs) {
    
        if (not ibf.interface) {
            continue;
        }

        for (const auto& itf : dg.interfaces) {
            if (*ibf.interface == itf) {
                std::complex<double> v{0.0, -ibf.length/(omega*MU0) };
                double sign = 1;
                if (sim.delta_gap_signs.size() > 0) {
                    sign = sim.delta_gap_signs[ibf.matrix_index];
                }
                context.V(ibf.matrix_index) += v*dg.voltage*sign;
            }
        }        
    }
}

void
update_rhs(simulation& sim, size_t ctx_number, const plane_wave& pw)
{
    auto& context = sim.contexts[ctx_number];
    auto& msh = sim.msh;

    auto propdir = pw.dir/pw.dir.norm();

    double freq = context.frequency;
    double omega = 2.0*M_PI*freq;
    auto kk = omega*std::sqrt(MU0*EPS0);

    for (const auto& ibf : sim.bfuncs) {

        std::complex<double> v{0.0, 0.0};

        const auto& iTminus = msh.triangles[ibf.itminus];
        const auto& iTplus = msh.triangles[ibf.itplus];

        auto wminus = 0.5*ibf.length/ibf.Aminus;
        auto wplus = 0.5*ibf.length/ibf.Aplus;

        auto mqps = integrate(msh, iTminus, sim.cfg.degree);
        for (auto& qp : mqps) {
            auto Rvec = qp.p - pw.srcpos;
            auto k = Rvec.to_eigen().dot(propdir)*kk;
            std::complex<double> jk{0.0, k};
            v -= wminus * qp.w * ibf.rho_minus(qp.p).to_eigen().dot(pw.E0)*std::exp(-jk);
        }

        auto pqps = integrate(msh, iTplus, sim.cfg.degree);
        for (auto& qp : pqps) {
            auto Rvec = qp.p - pw.srcpos;
            auto k = Rvec.to_eigen().dot(propdir)*kk;
            std::complex<double> jk{0.0, k};
            v += wplus * qp.w * ibf.rho_plus(qp.p).to_eigen().dot(pw.E0)*std::exp(-jk);
        }

        context.V(ibf.matrix_index) += v/std::complex<double>{0, -omega*MU0};        
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
        if (err.errtype == frico::meshing_error::multiple_triangles) {
            std::print(stderr,
"The mesh connectivity is not valid because an edge shared by more than two\n"
"triangles was detected. This does not permit to construct the RWG basis.\n"
"The tags of the involved surfaces are {}, {} and {}.\n", err.tminus_tag,
    err.tplus_tag, err.offending_tag);
            return false;
        }

        if (err.errtype == frico::meshing_error::lonely_edge) {
            std::print(stderr,
"Looks like that an edge from Curve {} is not on the boundary of any surface.\n"
"Probably you need to delete Curve {} from your GMSH geometry.\n",
    err.offending_tag, err.offending_tag);
            return false;
        }
    }

    std::println("Mesh information: ");
    std::println("        Vertices: {}", sim.msh.vertices.size());
    std::println("           Edges: {}", sim.msh.edges.size());
    std::println("           Cells: {}", sim.msh.triangles.size());
    std::println("           Beams: {}", sim.msh.beams.size());
    std::println("  Internal edges: {}", num_internal_edges(sim.msh));

    make_function_space(sim.msh, sim.bfuncs);

    return true;
}

bool init_sweep(simulation& sim, const frequency_range& freqs)
{
    size_t ctx_number = 0;
    for (double freq = freqs.start; freq <= freqs.end; freq += freqs.step) {
        freq_context context{};
        context.ctx_number = ctx_number++;
        context.frequency = freq;
        sim.contexts.push_back(context);
    }

    return true;
}

bool do_sweep(simulation& sim, const delta_gap& src)
{
    write_file_headers(sim, src);
    reorient_deltagap_edges(sim, src, sim.delta_gap_signs);

    for (size_t ctx_num = 0; ctx_num < sim.contexts.size(); ctx_num++) {
        run_context(sim, ctx_num, src);
        postpro_context(sim, ctx_num, src);

        /* Dump matrices, if needed*/
        if (sim.cfg.dump_matrices) {
            const auto& context = sim.contexts[ctx_num];
            std::string h5fn = "frico_" + std::to_string(ctx_num) + ".h5";
            H5Easy::File file(h5fn, H5Easy::File::Truncate);
            file.createDataSet("/frico/Z", context.Z);
            file.createDataSet("/frico/V", context.V);
        }

        /* dealloc matrix when we're done */
        sim.contexts[ctx_num].Z.resize(0,0);
    }
    return true;
}


void
postpro_context(const simulation& sim, size_t ctx_num, const plane_wave& dg)
{
    write_fields(sim, ctx_num);
}

bool do_sweep(simulation& sim, const plane_wave& src)
{
    for (size_t ctx_num = 0; ctx_num < sim.contexts.size(); ctx_num++) {
        run_context(sim, ctx_num, src);
        postpro_context(sim, ctx_num, src);

        /* Dump matrices, if needed*/
        if (sim.cfg.dump_matrices) {
            const auto& context = sim.contexts[ctx_num];
            std::string h5fn = "frico_" + std::to_string(ctx_num) + ".h5";
            H5Easy::File file(h5fn, H5Easy::File::Truncate);
            file.createDataSet("/frico/Z", context.Z);
            file.createDataSet("/frico/V", context.V);
        }

        /* dealloc matrix when we're done */
        sim.contexts[ctx_num].Z.resize(0,0);
    }
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

/*
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
*/






#if 0

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


    //make_radiation_diagrams(sim, ctx_number, {0,0,0}, 5.0);


    
    return true;
}

#endif






void write_file_headers(const simulation&, const plane_wave&)
{

}




void
postpro_context(simulation& sim, size_t ctx_number, const plane_wave& pw)
{
    std::ofstream ofs("rcs.txt", std::ios::out | std::ios::app);
    //auto [E, H] = eval_fields(sim, ctx_number, {0,0,-10});
    //ofs << sim.contexts[ctx_number].frequency << " " << 4*M_PI*100*(E.norm()*E.norm()) << std::endl;


    for (int i = 0; i < 359; i++) {
        double deg2rad = M_PI/180.0;
        double x = -10.0*std::sin(deg2rad*i);
        double z = -10.0*std::cos(deg2rad*i);
        auto [E, H] = eval_fields(sim, ctx_number, {x,0,z});

        auto xnorm = std::abs( E(0)*std::conj(E(0)) );
        auto ynorm = std::abs( E(1)*std::conj(E(1)) );
        auto znorm = std::abs( E(2)*std::conj(E(2)) );

        ezvec3 Ey{0.0, 1.0, 0.0};

        ofs << i << " " << 4*M_PI*100*( std::abs(E.dot(E)) )  << " ";
        ofs << 4*M_PI*100*xnorm << " ";
        ofs << 4*M_PI*100*ynorm << " ";
        ofs << 4*M_PI*100*znorm << " ";
        ofs << std::endl;
    }

    write_fields(sim, ctx_number);
}

}