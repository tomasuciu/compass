#ifndef TAUBIN_HPP
#define TAUBIN_HPP
#include "../fit.hpp"

namespace compass {
class TaubinSVD : AlgebraicFit<TaubinSVD> {
    public:
        using AlgebraicFit<TaubinSVD>::AlgebraicFit;

        TaubinSVD& fit (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data.rowwise() - mean;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            DesignMatrix mat = (DesignMatrix(centered.rows(), 3)
                    << centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            std::cout << mat << std::endl;
            return *this;
        }
};

class TaubinNewton : AlgebraicFit<TaubinNewton> {};

}
#endif /* TAUBIN_HPP */
