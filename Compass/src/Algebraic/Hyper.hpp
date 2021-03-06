#ifndef HYPER_HPP
#define HYPER_HPP

namespace compass {

class HyperSVD : public AlgebraicFit<HyperSVD> {
    friend class AlgebraicFit<HyperSVD>;

    public:
        HyperSVD() : AlgebraicFit<HyperSVD>() {}
        HyperSVD(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<HyperSVD>(data) {}

    protected:
        HyperSVD& compute (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = (Eigen::MatrixX<double>(data.rows(), data.cols())
                    << data.rowwise() - mean).finished();

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
                Eigen::RowVector4<double> R = mat.colwise().mean();
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

            AlgebraicFit::computeCircleParameters(solution);
            return *this;
        }

};

class HyperSimple : public AlgebraicFit<HyperSimple> {
    friend class AlgebraicFit<HyperSimple>;

    public:
        HyperSimple() : AlgebraicFit<HyperSimple>() {}
        HyperSimple(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<HyperSimple>(data) {}

    protected:
        HyperSimple& compute (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(centered.rows(), 4)
                    << Z, centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            Eigen::RowVectorX<double> S = mat.colwise().mean();
            Eigen::MatrixX<double> M = mat.transpose() * mat;

            Eigen::Matrix4<double> N = (Eigen::Matrix4<double>(4, 4)
                    << 8*S(0), 4*S(1), 4*S(2), 2, 4*S(1), 1, 0, 0, 4*S(2), 0, 1, 0, 2, 0, 0, 0).finished();

            Eigen::MatrixX<double> NM = N.inverse() * M;

            Eigen::EigenSolver<Eigen::MatrixXd> solver;
            solver.compute(NM);

            Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();
            Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

            // throw an exception and cancel execution
            if (eigenvalues(0) > 0) {
                std::cout << "Error, the smallest eigenvalue is positive" << std::endl;
            } else if (eigenvalues(1) < 0) {
                std::cout << "Error, the second smallest eigenvalue is negative" << std::endl;
            }

            Eigen::Vector4d AStar = eigenvectors.col(eigenvectors.cols() - 1);

            AlgebraicFit::computeCircleParameters(AStar, false);
            return *this;
        }
};
}

#endif /* HYPER_HPP */
