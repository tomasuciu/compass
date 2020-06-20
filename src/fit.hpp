#ifndef FIT_HPP
#define FIT_HPP

#include <numeric>
#include <algorithm>

#include "data.hpp"
#include "circle.hpp"

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

namespace compass {

struct Geometric {
static void LevenbergMarquardtFull() {}
static void LevenbergMarquardtReduced() {}
};

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>

struct Algebraic {

using ExtendedDesignMatrix = Eigen::Matrix<T, Eigen::Dynamic, 4>; // N x 4
using DesignMatrix = Eigen::Matrix<T, Eigen::Dynamic, 3>; // N x 3
using DataMatrix = Eigen::Matrix<T, 2, Eigen::Dynamic, Eigen::RowMajor>; // 2 X N
using DataMatrix3 = Eigen::Matrix<T, 3, Eigen::Dynamic>;

static void Kasa(const Eigen::Matrix<T, 2, Eigen::Dynamic, Eigen::RowMajor>& data) {
    auto mean = Eigen::Vector2<T>(data.row(0).mean(), data.row(1).mean());
    auto colwiseData = data.colwise();

    auto Mxy=0.0, Mxx=0.0, Myy=0.0, Mxz=0.0, Myz=0.0;
    std::for_each(colwiseData.begin(), colwiseData.end(), [&](const auto &column){

        auto Xi = column(0) - mean(0);
        auto Yi = column(1) - mean(1);
        auto Zi = std::pow(Xi, 2) + std::pow(Yi, 2);

        Mxx += std::pow(Xi, 2);
        Myy += std::pow(Yi, 2);
        Mxy += Xi * Yi;
        Mxz += Xi * Zi;
        Myz += Yi * Zi;
    });

    auto n = data.cols();
    Mxx /= n;
    Myy /= n;
    Mxy /= n;
    Mxz /= n;
    Myz /= n;

    Eigen::Matrix<T, 3, 3> augMatrix;
    augMatrix << Mxx, Mxy, 0,
                 Mxy, Myy, 0,
                 0,   0,   1;

    auto solVector = Eigen::Vector<T, 3>(-Mxz, -Myz, 0);

    Eigen::LDLT<Eigen::Matrix<T, 3, 3>> cholesky = augMatrix.ldlt();
    auto sol = cholesky.solve(solVector);
    auto B = -1 * sol(0)/2.0;
    auto C = -1 * sol(1)/2.0;

    Circle<T> fit;

    fit.a = B + mean(0);
    fit.b = C + mean(1);

    //TODO: Fix
    fit.radius = std::sqrt(std::pow(B, 2) + std::pow(C, 2));
    fit.sigma = 0; //TODO: Implement
    fit.i = 0;
    fit.j = 0;

    std::cout << fit << std::endl;
}

static void PrattNewton(const Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>& data) {
    auto mean = Eigen::Vector2<T>(data.row(0).mean(), data.row(1).mean());
    auto colwiseData = data.colwise();

    std::cout << mean << std::endl;
    auto Mxx=0.0, Myy=0.0, Mxy=0.0, Mxz=0.0, Myz=0.0, Mzz=0.0;

    std::for_each(colwiseData.begin(), colwiseData.end(), [&](const auto &column) {
        auto Xi = column(0) - mean(0);
        auto Yi = column(1) - mean(1);
        auto Zi = std::pow(Xi, 2) + std::pow(Yi, 2);

        Mxx += std::pow(Xi, 2);
        Myy += std::pow(Yi, 2);
        Mzz += std::pow(Zi, 2);
        Mxy += Xi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
    });

    int n = data.cols();

    Mxx /= n;
    Myy /= n;
    Mxy /= n;
    Mxz /= n;
    Myz /= n;
    Mzz /= n;

    // computing coefficients of characteristic polynomial

    double Mz = Mxx + Myy;
    double Cov_xy = Mxx*Myy - Mxy*Mxy;
    double Var_z = Mzz - Mz*Mz;

    double A2 = 4.0*Cov_xy - 3.0*Mz*Mz - Mzz;
    double A1 = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    double A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    double A22 = A2 + A2;

    double x, y, Dy, xnew, ynew;
    int iter, IterMAX=99;

    // finds root using newton's method
    for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)
        {
            Dy = A1 + x*(A22 + 16.*x*x);
            xnew = x - y/Dy;
            if ((xnew == x)||(!std::isfinite(xnew))) break;
            ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));
            if (abs(ynew)>=abs(y))  break;
            x = xnew;  y = ynew;
        }

    auto DET = x*x - x*Mz + Cov_xy;
    auto Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
    auto Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;

    std::cout << "X: " << Xcenter + mean(0) << std::endl;
    std::cout << "Y: " << Ycenter + mean(1) << std::endl;
    std::cout << "Radius: " << sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz + x + x) << std::endl;
}

static void PrattSVD(const Eigen::Matrix<T, 2, Eigen::Dynamic>& data) {

    Eigen::Vector2<T> mean{data.row(0).mean(), data.row(1).mean()};

    ExtendedDesignMatrix designMat(data.cols(), 4);

    auto designMatIt = designMat.rowwise().begin();
    for (const auto & col: data.colwise()) {
        auto Xi = col(0) - mean(0);
        auto Yi = col(1) - mean(1);
        auto Zi = std::pow(Xi, 2) + std::pow(Yi, 2);
        Eigen::Vector4<T> designMatRow{Zi, Xi, Yi, 1.0};

        *designMatIt = designMatRow;
        designMatIt = std::next(designMatIt);
    }

    Eigen::BDCSVD<Eigen::MatrixX<T>> svd(designMat, Eigen::ComputeThinU | Eigen::ComputeThinV);

    auto sigma = svd.singularValues();
    auto V = svd.matrixV();

    const double epsilon = 1.0e-13;
    if (sigma(3) < epsilon) {
        // singular case: the solution A is the right singular vector, i.e., fourth column of V
        auto A = V(3);
        // TODO: compute circle parameters.

        return;
    } else {
        // Y = V * sigma * V^T
        auto Y = V * sigma.asDiagonal() * V.transpose();

        Eigen::Matrix4<T> B;
        B << 0, 0, 0, -2,
             0, 1, 0,  0,
             0, 0, 1,  0,
            -2, 0, 0,  0;

        // compute Y * B^-1 * Y, where B is the constant constraint matrix
        auto result = Y * B.inverse() * Y;
        Eigen::EigenSolver<Eigen::MatrixX<T>> eigensolver;
        eigensolver.compute(result);

        Eigen::VectorX<T> eigen_values = eigensolver.eigenvalues().real();
        Eigen::MatrixX<T> eigen_vectors = eigensolver.eigenvectors().real();
        std::vector<std::pair<T, Eigen::VectorX<T>>> eigen_pairs;

        std::transform(eigen_values.begin(), eigen_values.end(), eigen_vectors.rowwise().begin(),
                std::back_inserter(eigen_pairs), [&](T val, Eigen::VectorX<T> vec) {
            return std::make_pair(val, vec);
        });

        // TODO: Instead of sorting, search for and return the smallest positive eigenpair
        // Consider using std::min_element()
        std::sort(std::begin(eigen_pairs), std::end(eigen_pairs),
                [&](const std::pair<T, Eigen::VectorX<T>>& a, const std::pair<T, Eigen::VectorX<T>>& b) {
            return std::get<0>(a) <= std::get<0>(b);
        });

        // Select the eigenpair (n, A_) with the smallest positive eigenvalue n

        int index = 0;
        std::pair<T, Eigen::VectorX<T>> eigenPairSmallestN;
        do {
            eigenPairSmallestN = eigen_pairs.at(index++);
        } while (std::get<0>(eigenPairSmallestN) <= 0.0);

        // compute A = Y^-1 A_
        auto A = Y.inverse() * std::get<1>(eigenPairSmallestN);

        // TODO: Why does this correction produce more stable results?
        auto A1_correction = -A(1)* 10;

        // compute circle parameters
        std::cout << "A correction : " << A1_correction << std::endl;
        auto a = -A(1) / A(0)/2.0 + mean(0);
        auto b = -A(2) / A(0)/2.0 + mean(1);
        auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) -4*A(0)*A(3)) / std::abs(A(0))/2.0;
        auto a_c = -A1_correction / A(0)/2.0 + mean(0);

        std::cout << "Center point : " << a << ", " << b << std::endl;
        std::cout << "Radius : " << radius << std::endl;
    }
}


static void Nievergelt(const DataMatrix& data) {
    using GolubDataMatrix = Eigen::Matrix<T, 3, Eigen::Dynamic>;
    GolubDataMatrix designMat(3, data.cols());

    Eigen::Vector2<T> mean{data.row(0).mean(), data.row(1).mean()};
    auto centered = data.colwise() - mean;

    auto centeredSquared = centered.cwiseProduct(centered);
    auto Mzz = centeredSquared.colwise().sum();

    auto Zmean = Mzz.mean();
    auto Zcentered = Mzz.array() - Zmean;

    designMat << Zcentered, centered;
    Eigen::BDCSVD<Eigen::MatrixX<T>> svd(designMat.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);

    auto V = svd.matrixV();
    auto A = V.col(2);
    auto A_4 = -Zmean * A(0);

    auto a = -A(1)/A(0)/2.0 + mean(0);
    auto b = -A(2)/A(0)/2.0 + mean(1);
    auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) - 4*A(0)*A_4)/std::abs(A(0))/2.0;
}


static void TaubinNewton(const DataMatrix& data) {

}

static void TaubinSVD(const DataMatrix& data) {
    DataMatrix3 matrix(3, data.cols());
    Eigen::Vector2<T> mean{data.row(0).mean(), data.row(1).mean()};
    auto centered = data.colwise() - mean;
    auto centeredSquared = centered.cwiseProduct(centered);
    auto Mzz = centeredSquared.colwise().sum();
    auto Zmean = Mzz.mean();
    auto Zi = Mzz.array() - Zmean;
    auto Zval = Zi/(2.0 * std::sqrt(Zmean));

    matrix << Zval, centered;

    Eigen::BDCSVD<Eigen::MatrixX<T>> svd(matrix.transpose(), Eigen::ComputeThinV);

    auto V = svd.matrixV();
    auto A = V.col(2);

    auto A0 = A(0)/(2.0*std::sqrt(Zmean));
    auto A_4 = -Zmean*A0;

    auto a = -A(1)/A0/2.0 + mean(0);
    auto b = -A(2)/A0/2.0 + mean(1);
    auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) - 4*A0 * A_4)/std::abs(A0)/2.0;

    std::cout << a << ", " << b << " | radius : " << radius << std::endl;
}

};

};
#endif
