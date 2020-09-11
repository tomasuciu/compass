#ifndef SPATH_HPP
#define SPATH_HPP

#include "../fit.hpp"
namespace compass {

template <class A>
class Spath : public GeometricFit<Spath<A>, A> {

    friend class GeometricFit<Spath<A>, A>;
    typedef GeometricFit<Spath<A>, A> Base;

    public:
        Spath<A>() : Base() {}
        Spath<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        Spath& compute (const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) {
            Eigen::MatrixXd centered = data.rowwise() - this->mean;
            std::cout << this->mean << std::endl;

            std::cout << initialGuess << std::endl;
            std::cout << centered << std::endl;

            return *this;
        }
};

}
#endif /* SPATH_HPP */
