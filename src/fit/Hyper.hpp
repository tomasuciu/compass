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

            Eigen::VectorXd solution;
            if (S.minCoeff() < 1e-12) {
                solution = V.col(3);
            } else {
                Eigen::RowVectorX<double> R = mat.colwise().mean();
                clamp(R);

                Eigen::MatrixX<double> W = V * S.asDiagonal() * V.transpose();

                Eigen::MatrixX<double> N = (Eigen::Matrix4<double>(4, 4)
                        << 8*R(0), 4*R(1), 4*R(2), 2, 4*R(1), 1, 0, 0, 4*R(2), 0, 1, 0, 2, 0, 0, 0).finished();

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

                solver.compute(W * N.inverse() * W);
                Eigen::MatrixXd eigenvectors = solver.eigenvectors();

                Eigen::Vector4d AStar = eigenvectors.col(1);
                Eigen::LDLT<Eigen::MatrixXd> cholesky = W.ldlt();
                solution = cholesky.solve(AStar);
            }

            auto a = -solution(1)/solution(0)/2.0 + mean(0);
            auto b = -solution(2)/solution(0)/2.0 + mean(1);

            auto radius = std::sqrt(std::pow(solution(1), 2) + std::pow(solution(2), 2) -4*solution(0)*solution(3)) / std::abs(solution(0))/2.0;
            std::cout << a << ", " << b << ", " << radius << std::endl;
            return *this;
        }

};
}

#endif /* HYPER_HPP */
