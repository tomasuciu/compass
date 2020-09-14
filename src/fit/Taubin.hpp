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
            circle.setParameters((-AR(1)/AR(0)/2) + mean(0), (-AR(2)/AR(0)/2) + mean(1), radius);

            return *this;
        }
};

class TaubinNewton : public AlgebraicFit<TaubinNewton> {
    friend class AlgebraicFit<TaubinNewton>;
    typedef AlgebraicFit<TaubinNewton> Base;

    public:
        TaubinNewton() : Base() {}
        TaubinNewton(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        TaubinNewton& compute (const Eigen::Ref<const DataMatrixD>& data) {
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

            double A3 = 4 * Mz;
            double A2 = -3 * Mz * Mz - Mzz;
            double A1 = Mzz * Mz + 4 * Cov_xy * Mz - NormedSymmetricMatrix(0, 1) * NormedSymmetricMatrix(0, 1)
                - NormedSymmetricMatrix(0, 2) * NormedSymmetricMatrix(0, 2) - Mz*Mz*Mz;

            double A0 = NormedSymmetricMatrix(0, 1) * NormedSymmetricMatrix(0, 1) * Myy + NormedSymmetricMatrix(0, 2) * NormedSymmetricMatrix(0, 2) * Mxx
                -Mzz * Cov_xy - 2* NormedSymmetricMatrix(0, 1) * NormedSymmetricMatrix(0, 2) * NormedSymmetricMatrix(1, 2) + Mz*Mz*Cov_xy;

            double A22 = A2 + A2;
            double A33 = A3 + A3 + A3;

            const int IterMax = 20;
            const double EPSILON = 1e-12;
            {
                int iter;
                double x , y;
                for (x = 0, y = 1e20, iter = 0; iter < IterMax; ++iter) {
                    double yold = y;
                    y = A0 + x*(A1 + x * (A2 + x*A3));
                    if (std::abs(y) > std::abs(yold)) {
                        // Newton - Taubin goes in the wrong direction!
                        x = 0;
                        break;
                    }
                    double Dy = A1 + x* (A22 + x*A33);
                    double xold = x;
                    x = xold - y/Dy;

                    if (std::abs((x - xold)/x) < EPSILON)
                        break;

                    if (x < 0) {
                        x = 0;
                    }
                }
                double DET = x*x - x*Mz + Cov_xy;
                double a = (NormedSymmetricMatrix(0, 1) * (Myy - x) - NormedSymmetricMatrix(0, 2) * NormedSymmetricMatrix(1, 2))/DET/2.0;
                double b = (NormedSymmetricMatrix(0, 2) * (Mxx - x) - NormedSymmetricMatrix(1, 0) * NormedSymmetricMatrix(1, 2))/DET/2.0;
                double rad = std::sqrt(a*a + b*b + Mz);
                std::cout << '(' << a + mean(0) << ", " << b + mean(1) << "), radius: " << rad << std::endl;
            }
            return *this;
        }
};

class TaubinNystromSVD : public AlgebraicFit<TaubinNystromSVD> {
    friend class AlgebraicFit<TaubinNystromSVD>;
    typedef AlgebraicFit<TaubinNystromSVD> Base;

    public:
        TaubinNystromSVD() : Base() {}
        TaubinNystromSVD(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        TaubinNystromSVD& compute (const Eigen::Ref<const DataMatrixD>& data) {
            // TODO: implement
            return *this;
        }
};

}
#endif /* TAUBIN_HPP */
