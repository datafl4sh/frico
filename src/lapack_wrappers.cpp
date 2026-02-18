#include <print>

#include <complex>
#include <vector>
#include <Eigen/Dense>

#include "lapack_wrappers.h"

extern "C" {
    void zsytrf_(const char* uplo, const int* n, std::complex<double>* a,
                 const int* lda, int* ipiv, std::complex<double>* work,
                 const int* lwork, int* info);

    void zsytrs_(const char* uplo, const int* n, const int* nrhs,
                 const std::complex<double>* a, const int* lda,
                 const int* ipiv, std::complex<double>* b,
                 const int* ldb, int* info);
}


bool solve(Eigen::MatrixXcd& A, const Eigen::VectorXcd& b,
    Eigen::VectorXcd& x) {
    assert(A.rows() == A.cols());
    assert(A.rows() == b.size());
    x = b;

    const int n = A.rows();
    const int nrhs = 1;
    int info;

    std::vector<int> ipiv(n);

    std::complex<double> workQuery;
    int lwork = -1;
    char uplo = 'U';

    zsytrf_(&uplo, &n, A.data(), &n, ipiv.data(),
            &workQuery, &lwork, &info);
    if(info != 0) {
        std::println(stderr, "ZSYTRF workspace query failed");
        return false;
    }

    lwork = static_cast<int>(std::real(workQuery));
    std::vector<std::complex<double>> work(lwork);

    zsytrf_(&uplo, &n, A.data(), &n, ipiv.data(),
            work.data(), &lwork, &info);
    if(info != 0) {
        std::println(stderr, "ZSYTRF factorization failed");
        return false;
    }
    
    zsytrs_(&uplo, &n, &nrhs, A.data(), &n, ipiv.data(),
            x.data(), &n, &info);

    if(info != 0) {
        std::println(stderr, "ZSYTRS solve failed");
        return false;
    }

    return true;
}







