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


static void Pratt(const Eigen::Ref<Eigen::Matrix2<T>>& data) {}

};

};
#endif
