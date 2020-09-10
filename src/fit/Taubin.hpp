#ifndef TAUBIN_HPP
#define TAUBIN_HPP
#include "../fit.hpp"
#include <Eigen/SVD>

namespace compass {
class TaubinSVD : public AlgebraicFit<TaubinSVD> {
    friend class AlgebraicFit<TaubinSVD>;

    public:
        TaubinSVD() : AlgebraicFit<TaubinSVD>() {}
        TaubinSVD(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<TaubinSVD>(data) {}

    protected:
        TaubinSVD& compute (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data.rowwise() - mean;
            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            Eigen::VectorX<double> Z0 = (Z.array() - Z.mean()) / (2 * std::sqrt(Z.mean()));

            DesignMatrix mat = (DesignMatrix(centered.rows(), 3)
                    << Z0, centered).finished();

            Eigen::BDCSVD<Eigen::MatrixX<double>> svd(mat, Eigen::ComputeThinV);
            Eigen::VectorX<double> A = svd.matrixV().col(2);
            A(0) = A(0) / (2 * std::sqrt(Z.mean()));

            Eigen::VectorX<double> AR = (Eigen::VectorX<double>(4) << A, -Z.mean() * A(0)).finished();

            // TODO: abstract away circle paramater calculation
            std::cout << (-AR(1)/AR(0)/2) + mean(0) << std::endl;
            std::cout << (-AR(2)/AR(0)/2) + mean(1) << std::endl;

            auto radius = std::sqrt(std::pow(AR(1), 2) + std::pow(AR(2), 2) - 4*AR(0) * AR(3))/std::abs(AR(0))/2.0;
            std::cout << radius << std::endl;

            return *this;
        }
};

class TaubinNewton : public AlgebraicFit<TaubinNewton> {
    public:
        using AlgebraicFit<TaubinNewton>::AlgebraicFit;

    protected:
        TaubinNewton& compute (const Eigen::Ref<const DataMatrixD>& data) {
            return *this;
        }
};

class TaubinNystromSVD : public AlgebraicFit<TaubinNystromSVD> {};

}
#endif /* TAUBIN_HPP */
