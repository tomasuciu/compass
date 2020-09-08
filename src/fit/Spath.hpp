#ifndef SPATH_HPP
#define SPATH_HPP
#include "../fit.hpp"
namespace compass {


    template <class A>
    class Spath : public GeometricFitExperimental<Spath<A>, A>{
        public:
            Spath<A>() : GeometricFitExperimental<Spath<A>, A>() {}
            Spath<A>(const Eigen::Ref<const DataMatrixD>& data) : GeometricFitExperimental<Spath<A>, A>(data) {}
    };
}
#endif /* SPATH_HPP */
