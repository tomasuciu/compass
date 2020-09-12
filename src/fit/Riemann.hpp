#ifndef RIEMANN_HPP
#define RIEMANN_HPP
#include "../fit.hpp"
namespace compass {

class Riemann : public SpecializedFit<Riemann> {
    friend class SpecializedFit<Riemann>;
    typedef SpecializedFit<Riemann> Base;

    public:
       Riemann() : Base() {}
       Riemann(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        Riemann& compute(const Eigen::Ref<const DataMatrixD>& data) {
            std::size_t rows = data.rows();
            Eigen::MatrixXd centered = data.rowwise() - this->mean;
            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            double factor = 2* std::sqrt(Z.sum()/rows);

            centered = centered.array() / factor;
            Z = Z.array() / factor / factor;

            Eigen::MatrixXd XYP = centered.array().colwise() / (Z.array() + 1);
            Eigen::VectorXd ZP = Z.array() / (Z.array() + 1);

            Eigen::VectorXd W = (Z.array() + 1) * (Z.array() + 1);
            Eigen::MatrixXd R = (Eigen::MatrixXd(XYP.rows(), 3)
                    << XYP, ZP).finished();

            Eigen::RowVector3d Rbar = (W.transpose() * R) / W.sum();
            Eigen::MatrixXd Rcentered = R.rowwise() - Rbar;
            Eigen::MatrixXd RScaled = Rcentered.array().colwise() * W.array();

            Eigen::MatrixXd scatter = RScaled.transpose() * Rcentered;

            Eigen::EigenSolver<Eigen::MatrixX<double>> solver;
            solver.compute(scatter, true);

            Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();
            Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

            if (eigenvalues.minCoeff() < 0) {
                std::cout << "Error, the smallest eigenvalue is negative." << std::endl;
                // throw new exception here
            }

            Eigen::RowVectorXd V = eigenvectors.col(eigenvectors.cols() - 1).transpose();
            Eigen::RowVector<double, 1> c = -Rbar * V.transpose();

            Eigen::RowVector2d Par = -V.head<2>().array() / 2 / (c(0) + V(2));
            Eigen::RowVector3d sol = (Eigen::RowVector3d()
                << Par, std::sqrt(Par*Par.transpose() - c(0)/(c(0) + V(2)))).finished();
            sol = sol.array() * factor;
            sol.head<2>() = sol.head<2>() + mean;

            std::cout << sol << std::endl;
            return *this;
        }
};

class RiemannAlgebraic : public SpecializedFit<RiemannAlgebraic>, public AlgebraicFit<RiemannAlgebraic> {
    friend class SpecializedFit<RiemannAlgebraic>;
    friend class AlgebraicFit<RiemannAlgebraic>;


};
}
#endif /* RIEMANN_HPP */
