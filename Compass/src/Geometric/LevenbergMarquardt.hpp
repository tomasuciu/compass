#ifndef LEVENBERG_MARQUARDT_HPP
#define LEVENBERG_MARQUARDT_HPP

#include "Compass/src/Core/fit.hpp"

namespace compass {
struct ObjectiveFunction {
    // typedef crazy tuple to be returned
    void computeIteration(double (*f)(Eigen::MatrixXd& data, const Circle<double>& current),
            const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& current) {

        Eigen::MatrixXd rescaled = data.rowwise() - current.getCenter();
        //auto [g, F, radius] = f(rescaled, current);

        std::cout << f(rescaled, current) << std::endl;
        // after shifted data matrix is returned, be sure to rescale by adding the mean colwise
    }
};


template<class A>
class LevenbergMarquardtFull : public GeometricFit<LevenbergMarquardtFull<A>, A> {
    friend class GeometricFit<LevenbergMarquardtFull<A>, A>;
    typedef GeometricFit<LevenbergMarquardtFull<A>, A> Base;

    public:
        LevenbergMarquardtFull(double lambda=1.0) : GeometricFit<LevenbergMarquardtFull<A>, A>(), lambda(lambda) {}
        LevenbergMarquardtFull(Eigen::Ref<DataMatrixD> data, Circle<double> guess, double lambda=1.0)
            : GeometricFit<LevenbergMarquardtFull<A>, A>(data, guess), lambda(lambda) {}

    protected:

        /*LevenbergMarquardtFull& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) {
            const double LAMBDA = 1.0;

            Circle<double> guess = initGuess;
            // intial computation of objective function and derivatives
            auto [J, g, F] = computeIteration(data, guess);

            const int iterMax = 50;
            const float epsilon = 0.000001;
            // main loop, each run is one iteration
            for (int i = 0; i < iterMax; ++i) {
                // computation of lambda
                while (true) {
                    Eigen::MatrixXd DelPar = (Eigen::MatrixXd(J.rows() + 3, J.cols())
                            << J, std::sqrt(lambda) * Eigen::MatrixXd::Identity(3, 3)).finished();

                    Eigen::Vector<double, 9> rhs = (Eigen::Vector<double, 9>()
                            << g, Eigen::Vector3d::Zero(3)).finished();

                    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr = DelPar.colPivHouseholderQr();

                    auto solution = qr.solve(rhs);
                    auto progress = computeNorm<double>(solution)/(computeNorm<double>(guess.getVector()) + epsilon);
                    std::cout << progress << std::endl;

                    exit(1);
                }
            }

            return *this;
        }*/

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

            Eigen::VectorX<double> g = D.array() - circle.getRadius();
            double F = std::pow(computeNorm<double>(g), 2);

            return std::make_tuple(J, g, F);
        }
};

template<class A>
class LevenbergMarquardtReduced : public GeometricFit<LevenbergMarquardtReduced<A>, A>, private ObjectiveFunction {
    friend class GeometricFit<LevenbergMarquardtReduced<A>, A>;
    typedef GeometricFit<LevenbergMarquardtReduced<A>, A> Base;

    public:
        LevenbergMarquardtReduced<A>() : Base() {}
        LevenbergMarquardtReduced<A>(const Eigen::Ref<const DataMatrixD>& data) {}

    protected:

        LevenbergMarquardtReduced<A>& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& initialGuess) {
            const double LAMBDA = 1.0;
            const double epsilon=0.000001;
            const int ITER_MAX = 50;

            double lambda_srt = std::sqrt(LAMBDA);

            computeIteration(&objective, data, initialGuess);

            for (int i = 0; i < ITER_MAX; ++i) {
                while (true) {

                    break;
                }

            }

            return *this;
        }

    private:

        //typedef std::tuple<Eigen::MatrixX<double>, Eigen::VectorX<double>
        // Note that the objective function modifies the data in place, reducing copies!
        static std::tuple<Eigen::VectorX<double>, double, double>
        objective(Eigen::MatrixXd& data, const Circle<double>& current) {

            Eigen::MatrixXd D = data.cwiseProduct(data);
            Eigen::VectorXd DSum = (D.rowwise().sum()).cwiseSqrt();

            data = data.array().colwise() / DSum.array();

            double radius = DSum.mean();
            Eigen::VectorXd g = DSum.array() - radius;
            double F = g.norm() * g.norm();

            return std::make_tuple(g, F, radius);
        }

};

}

#endif
