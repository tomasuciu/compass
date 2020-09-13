#ifndef INVERSION_HPP
#define INVERSION_HPP
#include "../fit.hpp"

namespace compass {

class InversionNonIterative : public SpecializedFitWithPole<InversionNonIterative> {
    friend class SpecializedFitWithPole<InversionNonIterative>;
    typedef SpecializedFitWithPole<InversionNonIterative> Base;

    public:
        InversionNonIterative() : Base() {}
        InversionNonIterative(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2<double>>& pole) : Base(data, pole) {}

    protected:
        InversionNonIterative& compute(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2<double>>& pole) {
            const size_t n = data.size();
            Eigen::MatrixX<double> shiftedToPole = data.rowwise() - pole.transpose();

            Eigen::VectorX<double> Z = (shiftedToPole.array().square()).rowwise().sum();
            DesignMatrix mat = (DesignMatrix(shiftedToPole.rows(), 3)
                    << Z, shiftedToPole).finished();

            Eigen::MatrixX<double> M = mat.transpose() * mat;

            Eigen::Matrix2<double> C = (Eigen::Matrix2<double>(2, 2)
                    << M.block<2, 2>(0, 0).determinant(), M({0, 1}, {0, 2}).determinant(), M({0, 1}, {0, 2}).determinant(), M({0, 2}, {0,2}).determinant()).finished();

            Eigen::EigenSolver<Eigen::MatrixXd> solver;
            solver.compute(C, true);

            Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();
            Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

            // Matrix must be positive semidefinite, i.e., eigenvalues must be nonnegative
            if (eigenvalues.real().minCoeff() < 0) {
                // throw an error and terminate execution, smallest eigenvalue is negative and thus not PSD
            }

            Eigen::Vector2d BC = eigenvectors.col(0);
            double A = -BC.transpose() * (M({1, 2}, {0}) / M(0, 0));

            Eigen::Vector2d center = (-BC / 2.0 / A) + pole;
            double radius = BC.norm() / 2.0 / std::abs(A);

            circle.setParameters(center(0), center(1), radius);
            return *this;
        }
};

class InversionIterative : public SpecializedFitWithPole<InversionIterative> {
    friend class SpecializedFitWithPole<InversionIterative>;
    typedef SpecializedFitWithPole<InversionIterative> Base;

    public:
        InversionIterative() : Base() {}
        InversionIterative(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2<double>>& pole) : Base(data, pole) {}

    protected:
        InversionIterative& compute(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2<double>>& pole) {
            const double EPSILON = 0.0000001;
            const int ITER_MAX = 50;

            Eigen::Vector2d pole1 = pole;
            Eigen::Vector2d pole2 = pole;

            Eigen::Vector3d par;
            for (int i = 0; i < ITER_MAX; ++i) {
                par.noalias() = InversionNonIterative(data, pole).getCircle().getVector();
                Eigen::Vector2d pole3 = 2 * par.head<2>() - pole2;

                if ((pole3 - pole1).norm() < EPSILON) {
                    break;
                }

                pole1 = pole2;
                pole2 = pole3;
            }
            // par contains the new parameters
            // TODO: initialize member circle accordingly
            return *this;
        }
};

}
#endif /* INVERSION_HPP */
