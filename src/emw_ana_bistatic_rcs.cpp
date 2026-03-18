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

#include "emw_postpro_common.h"
#include "emw_ana_bistatic_rcs.h"
#include "sources.h"
#include "utils.h"
namespace frico::maxwell {

bool
init_from_spec(const char *spec, bistatic_rcs_analysis& ana)
{
    auto toks = split(spec, ':');

    if ( toks.size() != 3 and toks.size() != 5 ) {
        std::println(stderr,
            "Invalid specification of bistatic RCS analysis (-b)");
        return false;
    }

    try {
        double deg2rad = M_PI/180.0;
        ana.radar_R = std::stod(toks[0]);
        ana.radar_theta = deg2rad*std::stod(toks[1]);
        ana.radar_phi = deg2rad*std::stod(toks[2]);

        if (toks.size() == 5) {
            ana.Etheta = std::stod(toks[4]);
            ana.Ephi = std::stod(toks[5]);
        }
    }
    catch (...)
    {
        std::println(stderr,
            "Error parsing specification of bistatic RCS analysis (-b)");
        return false;
    }

    return true;
}

static std::pair<double, double>
bistatic_RCS_unit_Ei(const ezvec3& Es, double theta_s, double phi_s, double R)
{
    ezvec3 theta_hat_s {
        std::cos(theta_s)*std::cos(phi_s),
        std::cos(theta_s)*std::sin(phi_s),
        -std::sin(theta_s) };

    ezvec3 phi_hat_s {
        -std::sin(phi_s),
        std::cos(phi_s),
        0.0 };

    std::complex<double> E_theta_s = theta_hat_s.dot(Es);
    std::complex<double> E_phi_s = phi_hat_s.dot(Es);

    double sigma_VV = 4.0 * M_PI * R * R * std::norm(E_theta_s);
    double sigma_HH = 4.0 * M_PI * R * R * std::norm(E_phi_s);

    return {sigma_VV, sigma_HH};
}

static void
postpro_context(simulation& sim, size_t ctx_num,
    const bistatic_rcs_analysis& ana, silo& db)
{
    if ( db.is_open() ) {
        std::string dirname = "sweep_step_" + std::to_string(ctx_num);
        auto old_dir = db.curdir().value();
        db.mkdir(dirname);
        db.chdir(dirname);
        write_fields(sim, ctx_num, db);
        db.chdir(old_dir);
    }

    std::string fname = "bistatic_rcs_" + std::to_string(ctx_num) + ".txt";

    std::ofstream ofsb(fname);
    std::println(ofsb, "# FRICO bistatic RCS computation");
    std::println(ofsb, "# Angle, VV_theta, HH_theta, VV_phi, HH_phi");
    for (int i = 0; i < 359; i++) {
        double deg2rad = M_PI/180.0;
        
        auto p_theta = sph2rect(ana.radar_R, deg2rad*i, 0.0);
        auto [Es_theta, Hs_theta] = eval_fields(sim, ctx_num, p_theta);
        auto [VV_theta, HH_theta] =
            bistatic_RCS_unit_Ei(Es_theta, deg2rad*i, 0.0, ana.radar_R);

        auto p_phi = sph2rect(ana.radar_R, 0.0, deg2rad*i);
        auto [Es_phi, Hs_phi] = eval_fields(sim, ctx_num, p_phi);
        auto [VV_phi, HH_phi] =
            bistatic_RCS_unit_Ei(Es_phi, 0.0, deg2rad*i, ana.radar_R);

        std::println(ofsb, "{} {} {} {} {}",
            i, VV_theta, HH_theta, VV_phi, HH_phi);
    }
}

bool
do_sweep(simulation& sim, const bistatic_rcs_analysis& ana)
{
    silo db;
    if (sim.cfg.silo_outfn) {
        db.open(sim.cfg.silo_outfn);
        db.mkdir("meshes");
        db.chdir("meshes");
        db.add_mesh("mesh", sim.msh);
    
        for (const auto& smp : sim.samplings.planes) {
            db.add_mesh(smp.name, smp.smpmsh);
        }
        db.chdir("/");
    }

    plane_wave pw(ana.radar_theta, ana.radar_phi, ana.Etheta, ana.Ephi);

    for (size_t ctx_num = 0; ctx_num < sim.contexts.size(); ctx_num++) {
        run_context(sim, ctx_num, pw);
        postpro_context(sim, ctx_num, ana, db);

        /* dealloc matrix when we're done */
        sim.contexts[ctx_num].Z.resize(0,0);
    }

    if ( db.is_open() ) {
        db.close();
    }

    return true;
}

} // namespace frico::maxwell