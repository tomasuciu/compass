#ifndef FIT_HPP
#define FIT_HPP

#include <algorithm>
#include "data.hpp"
#include "circle.hpp"
#include <Eigen/Dense>
#include <Eigen/Cholesky>

namespace compass {

struct Geometric {
static void LevenbergMarquardtFull() {}
static void LevenbergMarquardtReduced() {}
};

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>

struct Algebraic {

static void Kasa(const Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>& data) {
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

    auto n = data.size();
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


static void Pratt(const Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>& data) {
    auto mean = Eigen::Vector2<T>(data.row(0).mean(), data.row(1).mean());
    auto colwiseData = data.colwise();

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

    int n = data.size();
    Mxx /= n;
    Myy /= n;
    Mxy /= n;
    Mxz /= n;
    Myz /= n;
    Mzz /= n;

    // computing coefficients of characteristic polynomial
    auto Mz = Mxx + Myy;
    auto Cov_xy = Mxx*Myy - std::pow(Mxy, 2);
    auto Var_z = Mzz - std::pow(Mz, 2);

    auto A2 = 4.0*Cov_xy - 3.0*Mz*Mz - Mzz;
    auto A1 = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    auto A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    auto A22 = A2 + A2;

    int iter = 0;
    double x=0.0, y=0.0;
    for(x=0.0, y=0.0, iter=0; iter < 99; iter++) {

        double Dy = A1 + x*(A22 + 16.0 * std::pow(x, 2));
        double xnew = x - y/Dy;
        if ((xnew == x)||(!std::isfinite(xnew))) break;
        double ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));
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

};

};
#endif
