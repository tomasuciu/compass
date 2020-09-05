#ifndef PRATT_HPP
#define PRATT_HPP

#include <Eigen/SVD>
#include "../fit.hpp"

namespace compass {

class PrattSVD : public AlgebraicFit<PrattSVD> {
    public:
        PrattSVD() : AlgebraicFit<PrattSVD>() {}
        PrattSVD(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<PrattSVD>(data) {}

        PrattSVD& fit(const Eigen::Ref<const DataMatrixD>& data) {
            //this->mean = center<double>(data);
            Eigen::MatrixX<double> centered = data;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(centered.rows(), 4)
                    << Z, centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            Eigen::BDCSVD<Eigen::MatrixX<double>> svd{mat, Eigen::ComputeThinV};
            Eigen::MatrixX<double> V = svd.matrixV();
            Eigen::VectorX<double> S = svd.singularValues();

            if (S.minCoeff() < 1e-12) {
               // Eigen::VectorXd sol = V(
            } else {
                Eigen::MatrixX<double> W = V * S.asDiagonal();

                Eigen::MatrixX<double> Binv = (Eigen::Matrix4<double>(4, 4)
                        << 0, 0, 0, -0.5, 0, 1, 0, 0, 0, 0, 1, 0, -0.5, 0, 0, 0).finished();

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

                solver.compute(W.transpose() * Binv * W);
                Eigen::MatrixXd eigenvectors = solver.eigenvectors();
                Eigen::VectorXd eigenvalues = solver.eigenvalues();

                // (W^T * B^-1 * W) will only have one negative eigenvector, hence the smallest nonnegative
                // eigenvector must be the second (1st col)
                std::cout << eigenvalues << std::endl;
                Eigen::Vector4d smallestPositiveEigenvector = eigenvectors.col(1);

                Eigen::MatrixXd scaled = Eigen::Vector4d::Ones().cwiseQuotient(S).asDiagonal();

                Eigen::MatrixXd sol = V * scaled * smallestPositiveEigenvector;

                auto a = -sol(1)/sol(0)/2.0 + mean(0);
                auto b = -sol(2)/sol(0)/2.0 + mean(1);

                auto radius = std::sqrt(std::pow(sol(1), 2) + std::pow(sol(2), 2) -4*sol(0)*sol(3)) / std::abs(sol(0))/2.0;
                std::cout << a << ", " << b << ", " << radius << std::endl;
            }
            return *this;
        }
};

class PrattNewton : public AlgebraicFit<PrattNewton> {};

// Possibly using regula falsi with Illinois modification (i.e., introduces quasi-regularization parameters to
// scale function values according to iteration history)
class PrattRobust : public AlgebraicFit<PrattRobust> {};
}

#endif
