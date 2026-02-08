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

#include <iostream>
#include <fstream>
#include <print>

#include "eigen.h"


#include "output_silo.h"
#include "emw_solver.h"
#include "emw_postpro_common.h"
#include "emw_postpro_delta_gap.h"
namespace frico::maxwell {

/**
 * @brief Compute the radiation diagrams on the planes XY, YZ and XZ.
 * 
 * @param sim 
 * @param ctx_number 
 * @param pwr
 * @param gd
 * @return double
 */
void compute_radiation_diagrams(const simulation& sim, size_t ctx_number,
    double pwr, gain_data& gd)
{
    static const int steps = 360;

    gd.Gxy = ddvector::Zero(steps);
    gd.Gyz = ddvector::Zero(steps);
    gd.Gxz = ddvector::Zero(steps);

    #pragma omp parallel for
    for (int deg = 0; deg < steps; deg++) {
        double theta = deg*M_PI/180;
        double c = std::cos(theta);
        double s = std::sin(theta);
  
        /* XY */ {
            point Pxy{ gd.radius*c, gd.radius*s, 0.0 };
            auto R = norm(Pxy);
            auto [locE, locH] = eval_fields(sim, ctx_number, Pxy+gd.center);
            ezvec3 S = 0.5*locE.cross(locH.conjugate());
            double Prad = std::real(std::sqrt(S.dot(S)));
            double G = 4*M_PI*R*R*Prad/pwr;
            gd.Gxy(deg) = G;
        }

        /* YZ */ {
            point Pyz{ 0.0, gd.radius*c, gd.radius*s };
            auto R = norm(Pyz);
            auto [locE, locH] = eval_fields(sim, ctx_number, Pyz+gd.center);
            ezvec3 S = 0.5*locE.cross(locH.conjugate());
            double Prad = std::real(std::sqrt(S.dot(S)));
            double G = 4*M_PI*R*R*Prad/pwr;
            gd.Gyz(deg) = G;
        }

        /* XZ */ {
            point Pxz{ -gd.radius*c, 0.0, gd.radius*s };
            auto R = norm(Pxz);
            auto [locE, locH] = eval_fields(sim, ctx_number, Pxz+gd.center);
            ezvec3 S = 0.5*locE.cross(locH.conjugate());
            double Prad = std::real(std::sqrt(S.dot(S)));
            double G = 4*M_PI*R*R*Prad/pwr;
            gd.Gxz(deg) = G;
        }
    }

}

double
compute_max_gain(const gain_data& gd)
{
    auto maxGxy = gd.Gxy.maxCoeff();
    auto maxGyz = gd.Gyz.maxCoeff();
    auto maxGxz = gd.Gxz.maxCoeff();
    return std::max( {maxGxy, maxGyz, maxGxz} );
}

bool
write_radiation_diagrams(const std::string& filename, const gain_data& gd)
{
    std::ofstream ofs(filename);
    if ( not ofs.is_open() ) {
        std::println(stderr, "Can't open '{}' for writing", filename);
        return false;
    }
    for (int i = 0; i < 360; i++) {
        double magGxy = gd.Gxy(i);
        double magGyz = gd.Gyz(i);
        double magGxz = gd.Gxz(i);
        ofs << i << " " << magGxy << " " << 10*std::log10(magGxy)
                 << " " << magGyz << " " << 10*std::log10(magGyz)
                 << " " << magGxz << " " << 10*std::log10(magGxz)
                 << std::endl;
    }

    return true;
}

/**
 * @brief Compute current, impedance and power for a specific port modelled as
 *        a delta-gap
 * 
 * @param sim the simulation object
 * @param ctx_number context number
 * @param dg delta-gap source
 * @return port_values 
 */
port_values
compute_port_values(const simulation& sim,
    size_t ctx_number, const delta_gap& dg)
{
    auto& context = sim.contexts[ctx_number];

    port_values ret;

    for (const auto& ibf : sim.bfuncs) {
        if (not ibf.interface) {
            continue;
        }

        for (const auto& itf : dg.interfaces) {
            if (*ibf.interface == itf) {
                double sign = 1;
                if (sim.delta_gap_signs.size() > 0) {
                    sign = sim.delta_gap_signs[ibf.edge_index];
                }
                auto I = ibf.length*context.I(ibf.matrix_index)*sign;
                ret.I += I;
                ret.P += 0.5*dg.voltage*conj(I);
            }
        }   
    }
    
    ret.Z = dg.voltage/ret.I;

    return ret;
}

void
write_fields(const simulation& sim, size_t ctx_number)
{
    const auto& context = sim.contexts[ctx_number];

    std::string filename =
        sim.name + "_" + std::to_string(ctx_number) + ".silo";

    silo db;
    db.open(filename);
    db.add_mesh("mesh", sim.msh);
    db.add_variable("mesh", "J", context.tri_J, var_centering::zonal);
}

void write_file_headers(const simulation& sim, const delta_gap& dg)
{
    std::ofstream port_sweep_ofs("port_sweep.txt");
    std::println(port_sweep_ofs,
        "# Port sweep results, Z0 = {:.1f} Ohm",
        sim.cfg.Z0
    );
    std::println(port_sweep_ofs, "# step  frequency    Re(Z)    Im(Z)    SWR");
}

void
postpro_context(const simulation& sim, size_t ctx_num, const delta_gap& dg)
{
    const auto& context = sim.contexts[ctx_num];

    /* Computed port values */
    std::ofstream port_sweep_ofs("port_sweep.txt", std::ios_base::app);
    auto pv = compute_port_values(sim, ctx_num, dg);
    double swr = compute_swr(pv.Z, sim.cfg.Z0);
    std::println(port_sweep_ofs,
        "{:6}  {:10g}  {:10f}  {:10f}  {:10.3f}",
        ctx_num, context.frequency, std::real(pv.Z), std::imag(pv.Z), swr
    );
    port_sweep_ofs.close();

    std::println("Impedance Re/Im: ({:.4f},{:.4f}) Ohm",
        real(pv.Z), imag(pv.Z));

    std::println("Impedance abs/angle: {:.4f} Ohm, {:.4f} degrees",
        abs(pv.Z), 180*arg(pv.Z)/M_PI);
    
    std::println("SWR(Z0 = {:.1f} Ohm): {:.2f}",
        sim.cfg.Z0, swr);

    /* Radiation diagrams */
    std::string filename = "polar_" + std::to_string(ctx_num) + ".txt";
    gain_data rad_diags;
    compute_radiation_diagrams(sim, ctx_num, std::real(pv.P), rad_diags);
    write_radiation_diagrams(filename, rad_diags);

    auto maxG = compute_max_gain(rad_diags);
    std::println("Max gain: {:.2f} dB", 10*std::log10(maxG));
    
    /* Fields */
    write_fields(sim, ctx_num);
}

}