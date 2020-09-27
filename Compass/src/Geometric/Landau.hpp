#ifndef LANDAU_HPP
#define LANDAU_HPP

namespace compass {

template<class A>
class Landau : public GeometricFitImpl<Landau<A>, A> {
    friend GeometricFitImpl<Landau<A>, A>;
    typedef GeometricFitImpl<Landau<A>, A> Base;

    public:
        Landau<A>() : Base() {}
        Landau<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        Landau<A>& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) {
            Eigen::MatrixXd centered = data.rowwise() - this->mean;
            Eigen::RowVector3<double> ParNew = (Eigen::RowVector3<double>()
                    << (initialGuess.getCenter() - this->mean), 0).finished();

            const int MAX_ITER = 500;
            const double EPSILON = 0.00001;

            Eigen::RowVector3<double> ParOld;
            for (int i = 0; i < MAX_ITER; ++i) {
                ParOld.noalias() = ParNew;

                Eigen::MatrixXd Dxy = centered.rowwise() - ParOld.head<2>();
                Eigen::MatrixXd D = Dxy.cwiseProduct(Dxy);
                Eigen::VectorXd DSum = (D.rowwise().sum()).cwiseSqrt();

                Eigen::RowVector3d meanD = (Eigen::RowVector3d()
                        << (Dxy.array().colwise() / DSum.array()).colwise().mean(), 1).finished();

                meanD.head<2>() = meanD.head<2>().array() * -1;
                ParNew = meanD.array() * DSum.mean();

                double progress = (ParNew - ParOld).norm() / (ParOld.norm() + EPSILON);
                if (progress < EPSILON) {
                    break;
                }
            }
            ParOld.head<2>().noalias() = ParOld.head<2>() + this->mean;

            // ParOld is the solution vector [a, b, rad]
            std::cout << ParOld << std::endl;
            return *this;
        }
};
}
#endif /* LANDAU_HPP */
