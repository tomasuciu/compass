#ifndef SPATH_HPP
#define SPATH_HPP

namespace compass {

template <class A>
class Spath : public GeometricFitImpl<Spath<A>, A> {

    friend GeometricFitImpl<Spath<A>, A>;
    typedef GeometricFitImpl<Spath<A>, A> Base;

    public:
        Spath<A>() : Base() {}
        Spath<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        // TODO: Evaluate spurious copies to optimize space complexity
        Spath& compute (const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) {
            Eigen::MatrixXd centered = data.rowwise() - this->mean;
            Eigen::RowVector2<double> ParNew = initialGuess.getCenter() - this->mean;

            const int MAX_ITER = 500;
            Eigen::RowVector2<double> ParOld;

            for (int i = 0; i < MAX_ITER; ++i) {
                ParOld = ParNew;

                Eigen::MatrixXd Dxy = centered.rowwise() - ParOld;
                Eigen::MatrixXd D = Dxy.cwiseProduct(Dxy);
                Eigen::VectorXd DSum = (D.rowwise().sum()).cwiseSqrt();

                Eigen::MatrixXd rescaled = Dxy.array().colwise() / DSum.array();
                double Mu = rescaled.col(0).mean();
                double Mv = rescaled.col(1).mean();
                double Mr = DSum.mean();

                double radius = (Mu * ParOld(0) + Mv * ParOld(1) + Mr)/ (1 - std::pow(Mu, 2) - std::pow(Mv, 2));
                std::cout << radius << std::endl;

                ParNew(0) = -Mu * radius;
                ParNew(1) = -Mv * radius;

                // progress = (norm(ParNew - ParOld) / (norm(ParOld) + epsilon)
                // if (progress < epsilon) break;
                exit(1);
            }

            // recenter
            // Par = ParOld + mean;
            return *this;
        }
};

}
#endif /* SPATH_HPP */
