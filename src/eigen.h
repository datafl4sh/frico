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

#include <complex>

#ifdef HAVE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace frico {

using ddmatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using zdmatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
using ddvector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using zdvector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
using ddfield = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using zdfield = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 3>;

using dsparsematrix = Eigen::SparseMatrix<double>;
using zsparsematrix = Eigen::SparseMatrix<std::complex<double>>;

using edvec3 = Eigen::Matrix<double, 3, 1>;
using ezvec3 = Eigen::Matrix<std::complex<double>, 3, 1>;

} // namespace frico