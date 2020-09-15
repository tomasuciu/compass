#ifndef NIEVERGELT_HPP
#define NIEVERGELT_HPP

#include <eigen-master/Eigen/SVD>

#include "Compass/src/Core/fit.hpp"

namespace compass {

class Nievergelt : public AlgebraicFit<Nievergelt> {
    friend class AlgebraicFit<Nievergelt>;
    typedef AlgebraicFit<Nievergelt> Base;

    public:
        Nievergelt() : Base() {}
        Nievergelt(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}
    protected:
        Nievergelt& compute(const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data.rowwise() - mean;
            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            double Zmean = Z.mean();
            Eigen::VectorX<double> Zcentered = Z.array() - Z.mean();

            DesignMatrix mat = (DesignMatrix(centered.rows(), 3)
                    << Zcentered, centered).finished();
            std::cout << mat << std::endl;

            Eigen::BDCSVD<Eigen::MatrixX<double>> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

            auto V = svd.matrixV();
            auto A = V.col(2);
            auto A_4 = -Zmean * A(0);

            auto a = -A(1)/A(0)/2.0 + mean(0);
            auto b = -A(2)/A(0)/2.0 + mean(1);
            auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) - 4*A(0)*A_4)/std::abs(A(0))/2.0;
            std::cout << Circle<double>(a, b, radius) << std::endl;

            return *this;
        }
};
}

#endif /* NIEVERGELT_HPP */
