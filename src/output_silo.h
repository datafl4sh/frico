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

#pragma once

#include <string>

#include <Eigen/Dense>

#include "geom_mesh.h"
#include "silo.h"

namespace frico {

enum class var_centering {
    nodal,
    zonal
};

class silo {
    DBfile *db_;

public:
                silo();
    explicit    silo(const std::string&);
    explicit    silo(const char *);

                silo(const silo&) = delete;
    silo& operator=(const silo&) = delete;

    bool open(const std::string&);
    bool open(const char *);
    bool is_open() const;
    void close();

    bool add_mesh(const std::string& name, const mesh& msh);

    bool add_variable(const std::string& mesh_name,
        const std::string& var_name,
        const std::vector<double>& var,
        var_centering centering
    );

    bool add_variable(const std::string& mesh_name,
        const std::string& var_name,
        const Eigen::Matrix<double, Eigen::Dynamic, 1>& var,
        var_centering centering
    );

    bool add_variable(const std::string& mesh_name,
        const std::string& var_name,
        const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& var,
        var_centering centering
    );

    bool add_variable(const std::string& mesh_name,
        const std::string& var_name,
        const Eigen::Matrix<double, Eigen::Dynamic, 3>& var,
        var_centering centering
    );

    bool add_variable(const std::string& mesh_name,
        const std::string& var_name,
        const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 3>& var,
        var_centering centering
    );

    bool add_curve(const std::string& name, const std::vector<double>& x,
        const std::vector<double>& y
    );

    bool mkdir(const std::string& name);
    bool chdir(const std::string& name);
    std::optional<std::string> curdir(void) const;

    ~silo();
};

} // namespace frico
