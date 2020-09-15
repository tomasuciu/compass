#ifndef KUKUSH_MARKOVKSY_HUFFEL_HPP
#define KUKUSH_MARKOVKSY_HUFFEL_HPP

#include <eigen-master/Eigen/Eigenvalues>

#include "Compass/src/Core/fit.hpp"

namespace compass {
// TODO: create AlgebraicFitConsistent base class: derive KMH and KasaConsistent from it;
// Try and minimize dynamic matrix allocations
class KukushMarkovskyHuffel : public AlgebraicFit<KukushMarkovskyHuffel> {
    friend class AlgebraicFit<KukushMarkovskyHuffel>;
    typedef AlgebraicFit<KukushMarkovskyHuffel> Base;

    public:
        KukushMarkovskyHuffel() : Base() {}
        KukushMarkovskyHuffel(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        KukushMarkovskyHuffel& compute(const Eigen::Ref<const DataMatrixD>& data) {

        size_t rows = data.rows();

        // uncentered matrices
        Eigen::VectorX<double> Z = (data.array().square()).rowwise().sum();
        ExtendedDesignMatrix mat = (ExtendedDesignMatrix(rows, 4)
                << Z, data, Eigen::VectorXd::Ones(rows)).finished();

        DesignMatrix uncentered = (DesignMatrix(rows, 3)
                << data, Eigen::VectorXd::Ones(rows)).finished();

        Eigen::MatrixX<double> M0 = mat.transpose() * mat;

        Eigen::Matrix4<double> M1 = (Eigen::Matrix4<double>() <<
                8 * M0(0, 3), 4 * M0(1, 3), 4* M0(2, 3), 2 * rows, 4 * M0(1, 3), rows, 0, 0,
                4 * M0(2, 3), 0, rows, 0, 2 * rows, 0, 0, 0).finished();

        Eigen::Matrix4<double> M2 = (Eigen::Matrix4<double>() <<
                8 * rows, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();

        Eigen::MatrixX<double> centered = data.rowwise() - mean;
        Eigen::MatrixX<double> scatter = centered.transpose() * centered;

            Eigen::EigenSolver<Eigen::MatrixX<double>> solver;
            solver.compute(scatter, false);

            double Vmax = solver.eigenvalues().real().minCoeff();
            double Vmin = 0.0;

            float epsilon = 0.00001 * Vmax;
            Eigen::MatrixX<double> M;
            while (Vmax - Vmin > epsilon) {
                double V = (Vmin + Vmax)/2;
                M.noalias() = M0 - M1*V + V*V*M2;

                auto eig = M.eigenvalues();
                if (eig.real().minCoeff() > 0) {
                    Vmin = V;
                } else {
                    Vmax = V;
                }
            }

            solver.compute(M, true);
            Eigen::Index minIndex;
            double minEigenVal = solver.eigenvalues().real().minCoeff(&minIndex);
            Eigen::VectorXd minEigenVec = solver.eigenvectors().real().col(minIndex);

            auto a = -minEigenVec(1)/minEigenVec(0)/2.0;
            auto b = -minEigenVec(2)/minEigenVec(0)/2.0;
            auto rad = std::sqrt(a*a + b*b - minEigenVec(3)/minEigenVec(0));
            std::cout << Circle<double>(a, b, rad) << std::endl;
            return *this;
        }
};
}
#endif /*KUKUSH_MARKOVKSY_HUFFEL_HPP */
