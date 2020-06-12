#ifndef FIT_HPP
#define FIT_HPP

#include "data.hpp"
#include "circle.hpp"
#include <Eigen/Dense>

namespace compass {

struct Algebraic{};
struct Geometric{};

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>
static void Kasa(const Eigen::Ref<Eigen::Matrix2X<T>>& data) {
    double Xi, Yi, Zi;
    double Mxy, Mxx, Myy, Mxz, Myz;
    double B, C, G11, G12, G22, D1, D2;

    T meanX = data.row(0).mean();
    T meanY = data.row(1).mean();

    Mxx = Myy = Mxy = Mxz = Myz = 0;
    auto dataColWise = data.colwise();

    std::size_t index = 0;
    std::for_each(dataColWise.begin(), dataColWise.end(), [&](const auto &elem){
        // centered X and Y coordinates
        Xi = elem(0) - meanX;
        Yi = elem(1) - meanY;

        Zi = std::pow(Xi, 2) + std::pow(Yi, 2);

        // calculating moments
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

    // Naive implementation of Cholesky factorization; eventually use Eigen instead
    G11 = std::sqrt(Mxx);
    G12 = Mxy/G11;
    G22 = std::sqrt(Myy - std::pow(G12, 2));

    // calculating circle parameters
    D1 = Mxz/G11;
    D2 = (Myz - D1*G12)/G22;

    C = D2/G22/2.0;
    B = (D1 - G12*C)/G11/2.0;

    Circle<T> fit;

    fit.a = B + meanX;
    fit.b = C + meanY;
    fit.radius = std::sqrt(std::pow(B, 2) + std::pow(C, 2) + Mxx + Myy);
    fit.sigma = 0; //TODO: Implement utility to calculate RSS
    fit.i = 0;
    fit.j = 0;

    std::cout << fit << std::endl;
}

static void Pratt() {}
static void LevenbergMarquardtFull() {}
static void LevenbergMarquardtReduced() {}

}
#endif
