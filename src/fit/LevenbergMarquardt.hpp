#ifndef LEVENBERG_MARQUARDT_HPP
#define LEVENBERG_MARQUARDT_HPP

#include "../fit.hpp"

namespace compass {

class LevenbergMarquardtFull : public GeometricFit<LevenbergMarquardtFull> {
    public:
        LevenbergMarquardtFull(double lambda=1.0) : GeometricFit<LevenbergMarquardtFull>(), lambda(lambda) {}
        LevenbergMarquardtFull(Eigen::Ref<DataMatrixD> data, Circle<double> guess, double lambda=1.0)
            : GeometricFit<LevenbergMarquardtFull>(data, guess), lambda(lambda) {}

        LevenbergMarquardtFull& fit (Eigen::Ref<DataMatrixD> data, Circle<double> initGuess, double lambda=1.0) {
            Circle<double> guess = initGuess;

            // intial computation of objective function and derivatives
            auto [J, g, F] = computeIteration(data, guess);

            return *this;
        }

    private:
        double lambda;

        std::tuple<Eigen::MatrixX<double>, Eigen::VectorX<double>, double>
            computeIteration(Eigen::Ref<DataMatrixD> data, Circle<double> circle) {

            const size_t n = data.rows();
            Eigen::VectorX<double> Dx = data.col(0).array() - circle.getA();
            Eigen::VectorX<double> Dy = data.col(1).array() - circle.getB();

            Eigen::VectorX<double> D = (Dx.cwiseProduct(Dx) + Dy.cwiseProduct(Dy)).cwiseSqrt();
            Eigen::MatrixX<double> J = (Eigen::MatrixX<double>(n, 3)
                    << - Dx.cwiseQuotient(D), -Dy.cwiseQuotient(D), Eigen::VectorX<double>::Ones(n)).finished();

            std::cout << J << std::endl;

            Eigen::VectorX<double> g = D.array() - circle.getRadius();
            double F = std::pow(computeNorm<double>(g), 2);

            return std::make_tuple(J, g, F);
        }
};

class LevenbergMarquardtReduced : public GeometricFit<LevenbergMarquardtReduced> {};

}

#endif
