#ifndef NIEVERGELT_HPP
#define NIEVERGELT_HPP

namespace compass {

class Nievergelt : public AlgebraicFit<Nievergelt> {
    friend class AlgebraicFit<Nievergelt>;
    typedef AlgebraicFit<Nievergelt> Base;

    public:
        Nievergelt() : Base() {}
        Nievergelt(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        Nievergelt& compute(const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixXd centered = data.rowwise() - mean;
            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            double Zmean = Z.mean();

            Eigen::VectorX<double> Zcentered = Z.array() - Z.mean();
            DesignMatrix mat = (DesignMatrix(centered.rows(), 3)
                    << Zcentered, centered).finished();

            Eigen::BDCSVD<Eigen::MatrixX<double>> svd(mat, Eigen::ComputeThinV);

            Eigen::MatrixXd V = svd.matrixV();
            Eigen::Vector4d A = (Eigen::Vector4d()
                << V.col(V.cols() - 1), V.col(V.cols() - 1)(0) * -Zmean).finished();

            AlgebraicFit::computeCircleParameters(A);
            return *this;
        }
};
}

#endif /* NIEVERGELT_HPP */
