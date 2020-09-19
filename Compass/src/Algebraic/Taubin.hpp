#ifndef TAUBIN_HPP
#define TAUBIN_HPP

#include <random>

#include <eigen-master/Eigen/SVD>
#include <eigen-master/Eigen/QR>
#include <eigen-master/Eigen/Cholesky>

#include "Compass/src/Core/fit.hpp"

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

            AlgebraicFit::computeCircleParameters(AR);
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
                    y = A0 + x * (A1 + x * (A2 + x * A3));
                    if (std::abs(y) > std::abs(yold)) {
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
                computeCircleParameters(NormedSymmetricMatrix, x, y, Mxx, Myy, Mz, DET);
            }

            return *this;
        }

    private:
        void computeCircleParameters(const Eigen::Ref<const Eigen::MatrixXd>& NormedSymmetricMatrix,
                double x, double y, double Mxx, double Myy, double Mz, double DET) {

            double a = (NormedSymmetricMatrix(0, 1) * (Myy - x) - NormedSymmetricMatrix(0, 2) * NormedSymmetricMatrix(1, 2))/DET/2.0;
            double b = (NormedSymmetricMatrix(0, 2) * (Mxx - x) - NormedSymmetricMatrix(1, 0) * NormedSymmetricMatrix(1, 2))/DET/2.0;
            double rad = std::sqrt(a*a + b*b + Mz);
            circle.setParameters(a + mean(0), b + mean(1), rad);
        }

};

// TODO: Debug implementation!
class TaubinNystromSVD : public AlgebraicFit<TaubinNystromSVD> {
    friend class AlgebraicFit<TaubinNystromSVD>;
    typedef AlgebraicFit<TaubinNystromSVD> Base;

    public:
        TaubinNystromSVD() : Base() {}
        TaubinNystromSVD(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        // Note: This code is not configured to operate on Matrices of the DataMatrixD variety!
        TaubinNystromSVD& compute (const Eigen::Ref<const DataMatrixD>& data) {
            DataMatrix3 C(3, data.cols());

            auto centered = data.rowwise() - mean;
            auto centeredSquared = centered.cwiseProduct(centered);
            auto Mzz = centeredSquared.colwise().sum();
            auto Zmean = Mzz.mean();
            auto Zi = Mzz.array() - Zmean;
            auto Zval = Zi/(2.0 * std::sqrt(Zmean));

            C << Zval, centered;

            // generates a random gaussian matrix; abstract away eventually
            std::default_random_engine generator;
            std::normal_distribution<double> normal{};
            auto gaussian = [&] (double) {return normal(generator);};
            int l = 10;
            Eigen::MatrixXd omega = Eigen::MatrixXd::NullaryExpr(C.rows(), l, gaussian);

            Eigen::MatrixXd CTC = C * C.transpose();
            Eigen::MatrixXd X = CTC.colPivHouseholderQr().solve(omega);

            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
            qr.compute(X);

            Eigen::MatrixXd Q = qr.matrixQ();

            Eigen::MatrixXd Y = CTC.colPivHouseholderQr().solve(Q);

            Eigen::MatrixXd Z = Q.transpose() * Y;


            Eigen::LLT<Eigen::MatrixXd> cholesky = Z.llt();
            Eigen::MatrixXd G = cholesky.matrixL();
            Eigen::MatrixXd K = Y * G.inverse();

            Eigen::BDCSVD<Eigen::MatrixXd> svd(K, Eigen::ComputeFullU | Eigen::ComputeFullV);

            auto V = svd.matrixV();
            std::cout << "Nystrom SVD: \n" << V << std::endl;
            auto v = V.col(0);

            exit(1);
            auto A = V.col(2);

            auto A0 = A(0)/(2.0*std::sqrt(Zmean));
            auto A_4 = -Zmean*A0;

            auto a = -A(1)/A0/2.0 + mean(0);
            auto b = -A(2)/A0/2.0 + mean(1);
            auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) - 4*A0 * A_4)/std::abs(A0)/2.0;

            std::cout << a << ", " << b << " | radius : " << radius << std::endl;
            exit(1);
            return *this;
        }
};

}
#endif /* TAUBIN_HPP */
