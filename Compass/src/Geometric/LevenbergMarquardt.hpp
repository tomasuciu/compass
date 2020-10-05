#ifndef LEVENBERG_MARQUARDT_HPP
#define LEVENBERG_MARQUARDT_HPP

namespace compass {

template<class A>
class LevenbergMarquardtFull : public GeometricFitImpl<LevenbergMarquardtFull<A>, A> {
    friend GeometricFitImpl<LevenbergMarquardtFull<A>, A>;
    typedef GeometricFitImpl<LevenbergMarquardtFull<A>, A> Base;

    public:
        LevenbergMarquardtFull() : Base() {}
        LevenbergMarquardtFull(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        LevenbergMarquardtFull<A>& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& initialGuess) {
            const double LAMBDA = 1.0;

            const double epsilon=0.000001;
            const int ITER_MAX = 50;
            const int rows = data.rows();

            double lambda_srt = std::sqrt(LAMBDA);
            auto [J, g, F] = computeStep(data, initialGuess);
            Eigen::RowVector3d par = initialGuess.getVector().transpose();

            double progress;
            Eigen::RowVector3d ParTemp;
            Eigen::MatrixXd Jtemp;
            Eigen::VectorXd gtemp;
            double Ftemp;

            for (int i = 0; i < ITER_MAX; ++i) {
                double progress;
                while (true) {
                    Eigen::MatrixXd Del = (Eigen::MatrixXd(rows + 3, data.cols() + 1)
                            << J, Eigen::Matrix3d::Identity(3, 3).array() * lambda_srt).finished();

                    Eigen::VectorXd rhs = (Eigen::VectorXd(rows + 3)
                            << g, Eigen::Vector3d::Zero(3)).finished();

                    Eigen::Vector3d DelPar = Del.householderQr().solve(rhs);

                    progress = DelPar.norm() / (par.norm() + epsilon);
                    if (progress < epsilon) {
                        goto STOP;
                    }

                    ParTemp = par - DelPar.transpose();

                    std::tie(Jtemp, gtemp, Ftemp) = computeStep(data, ParTemp);

                    if (Ftemp < F && ParTemp(2) > 0) {
                        lambda_srt = lambda_srt / 2;
                        break;
                    } else {
                        lambda_srt = lambda_srt * 2;
                        continue;
                    }
                }

                if (progress < epsilon) {
                    break;
                }
                par = ParTemp;
                J = Jtemp;
                g = gtemp;
                F = Ftemp;
            }
            STOP:

            std::cout << par << std::endl;
            return *this;
        }

    private:
        //double lambda;

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

template<class A>
class LevenbergMarquardtReduced : public GeometricFitImpl<LevenbergMarquardtReduced<A>, A> {
    friend GeometricFitImpl<LevenbergMarquardtReduced<A>, A>;
    typedef GeometricFitImpl<LevenbergMarquardtReduced<A>, A> Base;

    public:
        LevenbergMarquardtReduced<A>() : Base() {}
        LevenbergMarquardtReduced<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        LevenbergMarquardtReduced<A>& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& initialGuess) {
            const double LAMBDA = 1.0;
            const double epsilon=0.000001;
            const int ITER_MAX = 50;

            const int rows = data.rows();

            double lambda_srt = std::sqrt(LAMBDA);

            Eigen::RowVector2d par = initialGuess.getCenter();
            auto [J, g, F, radius] = computeStep(data, par);

            for (int i = 0; i < ITER_MAX; ++i) {
                double progress;
                Eigen::RowVector2d ParTemp;
                Eigen::MatrixXd Jtemp;
                Eigen::VectorXd gtemp;
                double Ftemp, radiusTemp;

                while (true) {
                    Eigen::MatrixXd Del = (Eigen::MatrixXd(data.rows() + 2, data.cols())
                           << J, Eigen::Matrix2d::Identity(2, 2).array() * lambda_srt).finished();

                    Eigen::VectorXd rhs = (Eigen::VectorXd(data.rows() + 2)
                            << g, Eigen::Vector2d::Zero(2)).finished();

                    Eigen::VectorXd DelPar = Del.householderQr().solve(rhs);
                    ParTemp = par - DelPar.transpose();
                    progress = DelPar.norm() / (radius + par.norm() + epsilon);

                    if (progress < epsilon) {
                        goto STOP;
                    }

                    std::tie(Jtemp, gtemp, Ftemp, radiusTemp) = computeStep(data, ParTemp);

                    if (Ftemp < F) {
                        lambda_srt = lambda_srt / 2;
                        break;
                    } else {
                        lambda_srt = lambda_srt * 2;
                        continue;
                    }
                }

                par = ParTemp;
                J = Jtemp;
                g = gtemp;
                F = Ftemp;
                radius = radiusTemp;

                if (progress < epsilon)
                    break;
            }

            STOP:
            std::cout << par << std::endl;
            std::cout << radius << std::endl;

            return *this;
        }

    private:
        std::tuple<Eigen::MatrixXd, Eigen::VectorXd, double, double>
        computeStep(const Eigen::Ref<const DataMatrixD>& data, const Eigen::RowVector2d current) {
            Eigen::MatrixXd rescaled = data.rowwise() - current;

            Eigen::MatrixXd D = rescaled.cwiseProduct(rescaled);
            Eigen::VectorXd DSum = (D.rowwise().sum()).cwiseSqrt();

            Eigen::MatrixXd normalized = rescaled.array().colwise() / DSum.array();

            Eigen::RowVector2d normalizedMean = normalized.colwise().mean();

            double radius = DSum.mean();

            Eigen::MatrixXd J = -1 * normalized.array();
            J.noalias() = J.rowwise() + normalizedMean;

            Eigen::VectorXd g = DSum.array() - radius;
            double F = g.norm() * g.norm();
            return std::make_tuple(J, g, F, radius);
        }
};

} /* end compass */

#endif
