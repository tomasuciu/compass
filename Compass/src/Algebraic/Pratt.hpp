#ifndef PRATT_HPP
#define PRATT_HPP

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
            Eigen::Vector4d sol;

            if (S.minCoeff() < 1e-12) {
                sol = V.col(3);
            } else {
                Eigen::MatrixX<double> W = V * S.asDiagonal();

                Eigen::Matrix4<double> Binv = (Eigen::Matrix4<double>(4, 4)
                        << 0, 0, 0, -0.5, 0, 1, 0, 0, 0, 0, 1, 0, -0.5, 0, 0, 0).finished();

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

                solver.compute(W.transpose() * Binv * W);
                Eigen::MatrixXd eigenvectors = solver.eigenvectors();
                Eigen::Vector4d smallestPositiveEigenvector = eigenvectors.col(1);
                Eigen::MatrixXd scaled = Eigen::Vector4d::Ones().cwiseQuotient(S).asDiagonal();

                sol = V * scaled * smallestPositiveEigenvector;
            }

            AlgebraicFit::computeCircleParameters(sol);
            return *this;
        }
};

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
            double Cov_xy = Mxx * Myy - NormedSymmetricMatrix(1, 2) * NormedSymmetricMatrix(1, 2);

            double Var_z = Mzz - (Mz * Mz);

            double A2 = (4.0 * Cov_xy) - 3.0 * (Mz * Mz) - Mzz;

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

                double Xcenter = (NormedSymmetricMatrix(1,0) * (Myy - x) - NormedSymmetricMatrix(2, 0) * NormedSymmetricMatrix(1, 2))/DET/2.0;
                double Ycenter = (NormedSymmetricMatrix(2, 0) * (Mxx - x) - NormedSymmetricMatrix(1, 0) * NormedSymmetricMatrix(2, 1))/DET/2.0;
                double radius = sqrt(Xcenter * Xcenter + Ycenter * Ycenter + Mz + x + x);

                circle.setParameters(Xcenter + mean(0), Ycenter + mean(1), radius);
            }
            return *this;
        }
};

// Possibly using regula falsi with Illinois modification (i.e., introduces quasi-regularization parameters to
// scale function values according to iteration history)
class PrattRobust : public AlgebraicFit<PrattRobust> {};
}

#endif
