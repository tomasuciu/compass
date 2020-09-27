#ifndef GANDER_GOLUB_STREBEL_HPP
#define GANDER_GOLUB_STREBEL_HPP

namespace compass {
class GanderGolubStrebel : public AlgebraicFit<GanderGolubStrebel> {
    friend class AlgebraicFit<GanderGolubStrebel>;
    typedef AlgebraicFit<GanderGolubStrebel> Base;

    public:
        GanderGolubStrebel() : Base() {}
        GanderGolubStrebel(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        GanderGolubStrebel& compute (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::VectorX<double> Z = (Eigen::VectorX<double>(data.rows())
                    << (data.array().square()).rowwise().sum()).finished();

            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(data.rows(), 4)
                    << Z, data, Eigen::VectorXd::Ones(data.rows())).finished();

            Eigen::BDCSVD<Eigen::MatrixXd> SVD{mat, Eigen::ComputeThinV};
            Eigen::MatrixX<double> V = SVD.matrixV();

            AlgebraicFit::computeCircleParameters(SVD.matrixV().col(3));

            return *this;

        }
};
}
#endif /*GANDER_GOLUB_STREBEL_HPP */
