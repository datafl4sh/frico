#pragma once

#include <Eigen/Dense>

bool solve(Eigen::MatrixXcd& A, const Eigen::VectorXcd& b,
    Eigen::VectorXcd& x);