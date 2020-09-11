#ifndef PRATT_HPP
#define PRATT_HPP

#include <Eigen/SVD>
#include "../fit.hpp"

namespace compass {

class PrattSVD : public AlgebraicFit<PrattSVD> {
    friend class AlgebraicFit<PrattSVD>;

    public:
        PrattSVD() : AlgebraicFit<PrattSVD>() {}
        PrattSVD(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<PrattSVD>(data) {}

    protected:
        PrattSVD& compute(const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data.rowwise() - mean;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(centered.rows(), 4)
                    << Z, centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            Eigen::BDCSVD<Eigen::MatrixX<double>> svd{mat, Eigen::ComputeThinV};
            Eigen::MatrixX<double> V = svd.matrixV();
            Eigen::VectorX<double> S = svd.singularValues();

            if (S.minCoeff() < 1e-12) {
                // Solve idealized case here, TODO!
            } else {
                Eigen::MatrixX<double> W = V * S.asDiagonal();

                std::cout << W << std::endl;
                Eigen::MatrixX<double> Binv = (Eigen::Matrix4<double>(4, 4)
                        << 0, 0, 0, -0.5, 0, 1, 0, 0, 0, 0, 1, 0, -0.5, 0, 0, 0).finished();

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

                solver.compute(W.transpose() * Binv * W);
                Eigen::MatrixXd eigenvectors = solver.eigenvectors();

                Eigen::Vector4d smallestPositiveEigenvector = eigenvectors.col(1);
                std::cout << "\n\n" << smallestPositiveEigenvector << std::endl;

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
