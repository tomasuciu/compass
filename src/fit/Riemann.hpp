#ifndef RIEMANN_HPP
#define RIEMANN_HPP
#include "../fit.hpp"
namespace compass {

struct RiemannBase {
    void printFromBase() {
        std::cout << "printing from riemann base?" << std::endl;
    }
};

class Riemann : private RiemannBase, public SpecializedFit<Riemann> {
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

class RiemannAlgebraic : private RiemannBase, public SpecializedFit<RiemannAlgebraic> {
    friend class SpecializedFit<RiemannAlgebraic>;
    typedef SpecializedFit<RiemannAlgebraic> Base;

    public:
        RiemannAlgebraic() : Base () {}
        RiemannAlgebraic(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        RiemannAlgebraic& compute(const Eigen::Ref<const DataMatrixD>& data) {
            std::size_t rows = data.rows();
            Eigen::MatrixXd centered = data.rowwise() - this->mean;
            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();

            double factor = 2* std::sqrt(Z.sum()/rows);

            centered = centered.array() / factor;
            Z = Z.array() / factor / factor;

            Eigen::VectorX<double> Zp = Z.array() + 1;
            Eigen::VectorX<double> Zm = Z.array() - 1;

            double Zpp = Zp.transpose() * Zp;
            double Zpm = Zp.transpose() * Zm;
            double ZpX = Zp.transpose() * centered.col(0);
            double ZpY = Zp.transpose() * centered.col(1);

            Eigen::VectorX<double> A1 = (Zm.array() * Zpp - Zp.array()*Zpm) / 2.0;
            Eigen::VectorX<double> A2 = centered.col(0).array() * Zpp - Zp.array() * ZpX;
            Eigen::VectorX<double> A3 = centered.col(1).array() * Zpp - Zp.array() * ZpY;

            Eigen::MatrixX<double> AMat = (Eigen::MatrixX<double>(A1.rows(), 3)
                    << A1, A2, A3).finished();

            Eigen::BDCSVD<Eigen::MatrixXd> SVD{AMat, Eigen::ComputeThinV};
            Eigen::MatrixX<double> V = SVD.matrixV();

            Eigen::VectorX<double> P = V.col(V.cols() - 1);
            double Q = -2 * (P(0) * Zpm/2 + P(1) * ZpX + P(2) * ZpY) / Zpp;
            double A = (P(0) + Q) / 2.0;
            double D = (Q - P(0)) / 2.0;

            Eigen::RowVector3<double> solution = (Eigen::RowVector3<double>()
                    << -P(1)/A/2.0, -P(2)/A/2.0, std::sqrt((P.head<2>().transpose() * P.head<2>()) - 4*A*D)/2/std::abs(A)).finished();

            solution = solution.array() * factor;
            solution.head<2>() = solution.head<2>() + mean;

            std::cout << solution << std::endl;

            return *this;
        }
};
}
#endif /* RIEMANN_HPP */
