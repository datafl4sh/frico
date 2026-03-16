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

#include "emw_postpro_common.h"
#include "emw_postpro_rcs.h"

namespace frico::maxwell {

void write_file_headers(const simulation&, const plane_wave&)
{
    std::ofstream ofsm("monostatic_rcs.txt");
    ofsm << "# FRICO monostatic RCS computation" << std::endl;
    
    std::ofstream ofsb("bistatic_rcs.txt");
    
}

void
postpro_context(simulation& sim, size_t ctx_num, const plane_wave& pw)
{
    std::string dirname = "sweep_step_" + std::to_string(ctx_num);
    auto old_dir = sim.output_db.curdir().value();
    sim.output_db.mkdir(dirname);
    sim.output_db.chdir(dirname);
    write_fields(sim, ctx_num);
    sim.output_db.chdir(old_dir);

    std::ofstream ofsm("monostatic_rcs.txt", std::ios::out | std::ios::app);
    auto [E, H] = eval_fields(sim, ctx_num, {0,0,-10});
    ofsm << sim.contexts[ctx_num].frequency << " " << 4*M_PI*100*(E.norm()*E.norm()) << std::endl;

    std::string fname = "bistatic_rcs_" + std::to_string(ctx_num) + ".txt";

    std::ofstream ofsb(fname);
    ofsb << "# FRICO bistatic RCS computation" << std::endl;
    for (int i = 0; i < 359; i++) {
        double deg2rad = M_PI/180.0;
        double x = -10.0*std::sin(deg2rad*i);
        double z = -10.0*std::cos(deg2rad*i);
        auto [E, H] = eval_fields(sim, ctx_num, {x,0,z});

        auto rcs = 4*M_PI*100*( std::abs(E.dot(E)) );

        ofsb << i << " " <<  x << " " << z << " " << rcs << std::endl;
    }

}

} // namespace frico::maxwell