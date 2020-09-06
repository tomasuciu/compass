#ifndef HYPER_HPP
#define HYPER_HPP

#include "../fit.hpp"
#include "Eigen/SVD"

namespace compass {

class HyperSVD : public AlgebraicFit<HyperSVD> {
    public:
        HyperSVD() : AlgebraicFit<HyperSVD>() {}
        HyperSVD(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<HyperSVD>(data) {}

        HyperSVD& fit (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data.rowwise() - mean;
            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(centered.rows(), 4)
                    << Z, centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            Eigen::BDCSVD<Eigen::MatrixXd> SVD{mat, Eigen::ComputeThinV};
            Eigen::MatrixX<double> V = SVD.matrixV();
            Eigen::MatrixX<double> S = SVD.singularValues();

            if (S.minCoeff() < 1e-12) {
                // TODO Singular case, easy solution.
            } else {
                Eigen::RowVectorX<double> R = mat.colwise().mean();
                clamp(R);

                Eigen::MatrixX<double> W = V * S.asDiagonal();
                std::cout << R << std::endl;
            }

            return *this;
        }

};
}

#endif /* HYPER_HPP */
