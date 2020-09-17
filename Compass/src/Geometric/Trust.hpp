#ifndef TRUST_HPP
#define TRUST_HPP

#include "Compass/src/Core/fit.hpp"
#include <eigen-master/Eigen/QR>

namespace compass {
template <class A>
class Trust : public GeometricFit<Trust<A>, A> {
    friend class GeometricFit<Trust<A>, A>;
    typedef GeometricFit<Trust<A>, A> Base;

    public:
        Trust<A>() : Base() {}
        Trust<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        Trust& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& initialGuess) {
            double Delta = 1;
            const double EPSILON = 0.000001;
            int IterMax = 50;
            double c1 = 0.1, c2 = 0.25, c3 = 0.75;
            int nPar = 3;

            Eigen::RowVector3d par = initialGuess.getVector().transpose();
            auto [J, g, F] = computeStep(data, par);

            Eigen::VectorXd g0 = (Eigen::VectorXd(g.rows() + nPar)
                    << g, Eigen::VectorXd::Zero(nPar)).finished();

            for (int i = 0; i < IterMax; ++i) {
                while(true) {
                    double lambda = 0.0001;
                    Eigen::MatrixXd lhs = (Eigen::MatrixXd(J.rows() + nPar, J.cols())
                            << J, Eigen::Matrix3d::Identity().array() * sqrt(lambda)).finished();

                    Eigen::VectorXd h = lhs.householderQr().solve(g0);

                    if (h.norm() > (1 + c1) * Delta) {
                        double lambda_min = 0;
                        double lambda_max = 2;

                        while(true) {

                            lambda_max = lambda_max * lambda_max;
                            lhs = (Eigen::MatrixXd(J.rows() + nPar, J.cols())
                                    << J, Eigen::Matrix3d::Identity().array()  * sqrt(lambda)).finished();

                            h = lhs.householderQr().solve(g0);

                            if (h.norm() < Delta)
                                break;
                        }

                        lambda = 0.01;

                        while(true) {
                            // TODO
                            // find lambda given delta
                            lhs = (Eigen::MatrixXd(J.rows() + nPar, J.cols())
                                    << J, Eigen::Matrix3d::Identity().array()  * sqrt(lambda)).finished();
                            Eigen::HouseholderQR<Eigen::MatrixXd> qr(lhs.rows(), lhs.cols());
                            qr.compute(lhs);

                            Eigen::MatrixXd q = qr.matrixQR();
                        }


                    }
                }
            }

            return *this;
        }
    private:

        std::tuple<Eigen::MatrixX<double>, Eigen::VectorX<double>, double>
        computeStep(const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& current) {
            const size_t rows = data.rows();

            Eigen::MatrixXd rescaled = data.rowwise() - current.getCenter();

            Eigen::MatrixXd D = rescaled.cwiseProduct(rescaled);
            Eigen::VectorXd DSum = (D.rowwise().sum()).cwiseSqrt();

            Eigen::MatrixXd normalized = rescaled.array().colwise() / DSum.array();


            Eigen::MatrixXd J = (Eigen::MatrixXd(rows, data.cols() + 1)
                    << normalized.array(), Eigen::VectorXd::Ones(rows)).finished();
            J = J.array() * -1;
            //std::cout << J << std::endl;
            Eigen::VectorXd g = DSum.array() - current.getRadius();
            double F = g.norm() * g.norm();

            return std::make_tuple(J, g, F);
        }

        std::tuple<Eigen::MatrixX<double>, Eigen::VectorX<double>, double>
        computeStep(const Eigen::Ref<const DataMatrixD>& data, const Eigen::RowVector3d current) {
            Circle<double> temp(current(0), current(1), current(2));
            return computeStep(data, temp);
        }
};

}

#endif /* TRUST_HPP */
