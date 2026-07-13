/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025,2026 Matteo Cicuttin - IV3IWE
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

#include "emw_ana_antenna.h"
#include "emw_postpro_common.h"

namespace frico::maxwell {

enum class gain_plane {
    xy,
    yz,
    xz
};

static void compute_gain(const simulation& sim, size_t ctx_number,
    double pwr, double R, gain_plane plane, std::vector<double>& Gv)
{
    static const int steps = 360;
    Gv.resize(steps);

    #pragma omp parallel for
    for (int deg = 0; deg < steps; deg++) {
        double theta = deg*M_PI/180;
        double c = std::cos(theta);
        double s = std::sin(theta);

        point P;
        switch (plane) {
            case gain_plane::xy:
                P = point{ R*c, R*s, 0.0 };
                break;
            case gain_plane::yz:
                P = point{ 0.0, R*c, R*s };
                break;
            case gain_plane::xz:
                P = point{ -R*c, 0.0, R*s };
                break;
        }
  
        auto R = norm(P);
        auto [locE, locH] = eval_fields(sim, ctx_number, P);
        ezvec3 S = 0.5*locE.cross(locH.conjugate());
        double Prad = std::real(std::sqrt(S.dot(S)));
        double G = 4*M_PI*R*R*Prad/pwr;
        Gv[deg] = G;
    }
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
static void
compute_port_values(const simulation& sim,
    size_t ctx_number, antenna_analysis& anta)
{
    auto& context = sim.contexts[ctx_number];
    auto& antdata = anta.antdata[ctx_number];

    for (const auto& ibf : sim.bfuncs) {
        if (not ibf.interface) {
            continue;
        }

        for (const auto& itf : anta.dgap.interfaces) {
            if (*ibf.interface == itf) {
                double sign = 1;
                if (sim.delta_gap_signs.size() > 0) {
                    sign = sim.delta_gap_signs[ibf.matrix_index];
                }
                auto I = ibf.length*context.I(ibf.matrix_index)*sign;
                antdata.port.current += I;
                antdata.port.power += 0.5*anta.dgap.voltage*conj(I);
            }
        }   
    }
    
    antdata.port.voltage = anta.dgap.voltage;
    antdata.port.impedance = anta.dgap.voltage/antdata.port.current;
    antdata.port.gamma = gamma(antdata.port.impedance, sim.cfg.Z0);
    antdata.port.swr = swr(antdata.port.gamma);
}

static void
postpro_context(simulation& sim, size_t ctx_num,
    antenna_analysis& anta, silo& db)
{
    auto& antdata = anta.antdata[ctx_num];

    /* Computed port values */
    compute_port_values(sim, ctx_num, anta);

    std::println("Impedance Re/Im: ({:.4f},{:.4f}) Ohm",
        real(antdata.port.impedance), imag(antdata.port.impedance));

    std::println("Impedance abs/angle: {:.4f} Ohm, {:.4f} degrees",
        abs(antdata.port.impedance), 180*arg(antdata.port.impedance)/M_PI);
    
    std::println("SWR(Z0 = {:.1f} Ohm): {:.2f}",
        sim.cfg.Z0, antdata.port.swr);

    /* Radiation diagrams */
    auto pwr = real(antdata.port.power);
    compute_gain(sim, ctx_num, pwr, anta.rdiag_dist,
        gain_plane::xy, antdata.gain_XY);
    compute_gain(sim, ctx_num, pwr, anta.rdiag_dist,
        gain_plane::yz, antdata.gain_YZ);
    compute_gain(sim, ctx_num, pwr, anta.rdiag_dist,
        gain_plane::xz, antdata.gain_XZ);

    auto maxG_XY = *std::max_element(antdata.gain_XY.begin(), antdata.gain_XY.end());
    auto maxG_YZ = *std::max_element(antdata.gain_YZ.begin(), antdata.gain_YZ.end());
    auto maxG_XZ = *std::max_element(antdata.gain_XZ.begin(), antdata.gain_XZ.end());
    auto maxG = std::max({ maxG_XY, maxG_YZ, maxG_XZ });
    std::println("Max gain: {:.2f} dB", 10*std::log10(maxG));
    
    /* Fields */
    if ( db.is_open() ) {
        std::string dirname = "sweep_step_" + std::to_string(ctx_num);
        auto old_dir = db.curdir().value();
        db.mkdir(dirname);
        db.chdir(dirname);
        write_fields(sim, ctx_num, db);
        db.chdir(old_dir);
    }
}

static bool
write_sweep_data(const std::string& filename, const simulation& sim,
    const antenna_analysis& anta)
{
    std::ofstream port_sweep_ofs(filename);
    std::println(port_sweep_ofs,
        "# Port sweep results, Z0 = {:.1f} Ohm",
        sim.cfg.Z0
    );
    std::println(port_sweep_ofs, "# step  frequency    Re(Z)    Im(Z)    SWR");

    auto nctxs = sim.contexts.size();
    for (size_t ctx_num = 0; ctx_num < nctxs; ctx_num++) {
        const auto& context = sim.contexts[ctx_num];
        const auto& antdata = anta.antdata[ctx_num];
        std::println(port_sweep_ofs, "{} {} {} {} {}", ctx_num,
            context.frequency, antdata.port.impedance.real(),
            antdata.port.impedance.imag(), antdata.port.swr);
    }

    return true;
}

static bool
write_radiation_diagrams(const std::string& filename, const simulation& sim,
    size_t ctx_num, const antenna_analysis& anta)
{
    std::ofstream ofs(filename);
    if ( not ofs.is_open() ) {
        std::println(stderr, "Can't open '{}' for writing", filename);
        return false;
    }

    const auto& antdata = anta.antdata[ctx_num];

    std::println(ofs, "# FRICO radiation patterns. Frequency = {}",
        sim.contexts[ctx_num].frequency);

    std::println(ofs, "# angle XY_lin XY_dB YZ_lin YZ_dB XZ_lin XZ_dB");

    for (int i = 0; i < 360; i++) {
        double magGxy = antdata.gain_XY[i];
        double magGyz = antdata.gain_YZ[i];
        double magGxz = antdata.gain_XZ[i];

        std::println(ofs, "{} {} {} {} {} {} {}", i,
            magGxy, 10.0*std::log10(magGxy),
            magGyz, 10.0*std::log10(magGyz),
            magGxz, 10.0*std::log10(magGxz)
        );
    }

    return true;
}


static bool
write_hdf5(const std::string& filename, const simulation& sim,
    const antenna_analysis& anta)
{
    auto nctxs = sim.contexts.size();
    const std::string basepath = "/frico/antenna_analysis/";
    H5Easy::File file(filename, H5Easy::File::Truncate);
    file.createDataSet(basepath+"rdiag_dist", anta.rdiag_dist);
    file.createDataSet(basepath + "sweep_steps", nctxs);
    file.createDataSet(basepath + "Z0", sim.cfg.Z0);
    for (size_t ctx_num = 0; ctx_num < nctxs; ctx_num++) {
        const auto& context = sim.contexts[ctx_num];
        const auto& antdata = anta.antdata[ctx_num];
        const std::string path = basepath + std::to_string(ctx_num) + "/";
        file.createDataSet(path + "frequency", context.frequency);
        file.createDataSet(path + "gain_XY", antdata.gain_XY);
        file.createDataSet(path + "gain_YZ", antdata.gain_YZ);
        file.createDataSet(path + "gain_XZ", antdata.gain_XZ);
        file.createDataSet(path + "port/voltage", antdata.port.voltage);
        file.createDataSet(path + "port/current", antdata.port.current);
        file.createDataSet(path + "port/power", antdata.port.power);
        file.createDataSet(path + "port/impedance", antdata.port.impedance);
        file.createDataSet(path + "port/gamma", antdata.port.gamma);
        file.createDataSet(path + "port/swr", antdata.port.swr);
    }

    return true;
}

bool do_sweep(simulation& sim, antenna_analysis& anta)
{
    reorient_deltagap_edges(sim, anta.dgap, sim.delta_gap_signs);

    auto nctxs = sim.contexts.size();
    anta.antdata.resize(nctxs);

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

    for (size_t ctx_num = 0; ctx_num < nctxs; ctx_num++) {
        run_context(sim, ctx_num, anta.dgap);
        postpro_context(sim, ctx_num, anta, db);

        /* Dump matrices, if needed*/
        if (sim.cfg.dump_matrices) {
            const auto& context = sim.contexts[ctx_num];
            std::string h5fn = "frico_linear_system_" + std::to_string(ctx_num) + ".h5";
            H5Easy::File file(h5fn, H5Easy::File::Truncate);
            file.createDataSet("/frico/Z", context.Z);
            file.createDataSet("/frico/V", context.V);
            file.createDataSet("/frico/I", context.I);
        }

        /* dealloc matrix when we're done */
        sim.contexts[ctx_num].Z.resize(0,0);
    }

    if (sim.cfg.write_text_outfiles) {
        /* Write text files with radiation diagrams */
        for (size_t ctx_num = 0; ctx_num < nctxs; ctx_num++) {
            std::string fname = "polar_" + std::to_string(ctx_num) + ".txt";
            write_radiation_diagrams(fname, sim, ctx_num, anta);
        }

        /* Write text file with computed port parameters */
        write_sweep_data("port_sweep.txt", sim, anta);
    }

    if (sim.cfg.h5_outfn) {
        write_hdf5(sim.cfg.h5_outfn, sim, anta);
    }


    if ( db.is_open() ) {
        db.close();
    }

    return true;
}

} // namespace frico::maxwell