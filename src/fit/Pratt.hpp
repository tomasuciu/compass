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

                //std::cout << W << std::endl;
                Eigen::MatrixX<double> Binv = (Eigen::Matrix4<double>(4, 4)
                        << 0, 0, 0, -0.5, 0, 1, 0, 0, 0, 0, 1, 0, -0.5, 0, 0, 0).finished();

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

                solver.compute(W.transpose() * Binv * W);
                Eigen::MatrixXd eigenvectors = solver.eigenvectors();

                Eigen::Vector4d smallestPositiveEigenvector = eigenvectors.col(1);
                //std::cout << "\n\n" << smallestPositiveEigenvector << std::endl;

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

// TODO:
// 1 - Replace dynamic matrices with constant-size matrices for improved performance
// 2 - Replace instances of std::pow(X, 2) with X * X
// 3 - Write method to compute circle parameters, linking with .getCircle()
// 4 - Design and implement PrattBase mixin for creating design matrices
class PrattNewton : public AlgebraicFit<PrattNewton> {
    friend class AlgebraicFit<PrattNewton>;
    typedef AlgebraicFit<PrattNewton> Base;

    public:
        PrattNewton() : Base() {}
        PrattNewton(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        PrattNewton& compute (const Eigen::Ref<const DataMatrixD>& data) {
            const std::size_t N = data.rows();
            Eigen::MatrixX<double> centered = data.rowwise() - mean;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(N, 4)
                    << Z, centered, Eigen::VectorXd::Ones(N)).finished();

            Eigen::MatrixX<double> NormedSymmetricMatrix = (mat.transpose() * mat).array() / N;
            clamp(NormedSymmetricMatrix);

            Eigen::Vector4<double> NormedDiagonal = NormedSymmetricMatrix.diagonal();

            double Mzz = NormedDiagonal(0), Mxx = NormedDiagonal(1), Myy = NormedDiagonal(2);
            double Mz = Mxx + Myy;
            double Cov_xy = Mxx * Myy - std::pow(NormedSymmetricMatrix(1, 2), 2);

            double Var_z = Mzz - std::pow(Mz, 2);

            double A2 = (4.0 * Cov_xy) - 3.0 * std::pow(Mz, 2) - Mzz;

            double A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - std::pow(NormedSymmetricMatrix(0, 1), 2) - std::pow(NormedSymmetricMatrix(0, 2), 2);
            double A0 = NormedSymmetricMatrix(0, 1) * (NormedSymmetricMatrix(0, 1)*Myy - NormedSymmetricMatrix(0, 2)*NormedSymmetricMatrix(2, 1))
                + NormedSymmetricMatrix(0, 2)*(NormedSymmetricMatrix(0, 2)*Mxx - NormedSymmetricMatrix(1, 0)*NormedSymmetricMatrix(2, 1)) - Var_z*Cov_xy;

            const int IterMAX = 99;

            {
                int iter;
                double x = 0, y = A0;
                for (x = 0, y = A0, iter = 0; iter < IterMAX; ++iter) {
                    double Dy = A1 + x * (2 * A2 + 16.0 * (x * x));
                    double xnew = x - y/Dy;

                    if ((xnew == x) || (!std::isfinite(xnew)))
                        break;

                    double ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));

                    if (abs(ynew) >= abs(y))
                        break;
                    x = xnew;  y = ynew;
                }

            double DET = x*x - x*Mz + Cov_xy;

            double Xcenter = (NormedSymmetricMatrix(1,0)*(Myy - x) - NormedSymmetricMatrix(2, 0)*NormedSymmetricMatrix(1, 2))/DET/2.0;
            double Ycenter = (NormedSymmetricMatrix(2, 0)*(Mxx - x) - NormedSymmetricMatrix(1, 0)*NormedSymmetricMatrix(2, 1))/DET/2.0;

            std::cout << "X: " << Xcenter + mean(0) << std::endl;
            std::cout << "Y: " << Ycenter + mean(1) << std::endl;
            std::cout << "Radius: " << sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz + x + x) << std::endl;
            }
            return *this;
        }
};

// Possibly using regula falsi with Illinois modification (i.e., introduces quasi-regularization parameters to
// scale function values according to iteration history)
class PrattRobust : public AlgebraicFit<PrattRobust> {};
}

#endif
