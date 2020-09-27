#ifndef FORWARD_DECLARATIONS_HPP
#define FORWARD_DECLARATIONS_HPP

namespace compass {

//TODO: Add macro to check for eigen and the include the necessary header files

template<typename Derived> struct FitCRTP;
template<typename Derived> class FitBase;

template<typename Derived> class AlgebraicFit;

class TaubinNewton;
class TaubinSVD;
class TaubinNystromSVD;
class PrattNewton;
class PrattSVD;
class PrattRobust;
class Nievergelt;
class KukushMarkovskyHuffel;
class Kasa;
class KasaConsistent;
class HyperSVD;
class HyperSimple;
class GanderGolubStrebel;

template<class...>
struct False : std::integral_constant<bool, false>{};

template<typename Derived, class A, bool Enable>
class GeometricFit;

template<typename Derived, class A>
class GeometricFit<Derived, A, true>;

template<typename Derived, class A>
class GeometricFit<Derived, A, false>;

template<typename Derived, class A>
using GeometricFitImpl = GeometricFit<Derived, A, (std::is_base_of<AlgebraicFit<A>, A>::value)>;


template<typename Derived> class GeometricFitWithGuess;
template<typename Derived> class GeometricFitGuessIndiscriminate;


template<typename Derived> class SpecializedFit;
template<typename Derived> class SpecializedFitWithPole;
template<typename Derived> class SpecializedFitRandomPole;

// test geometric methods
template<class A> class Landau;
template<class A> class Trust;
template<class A> class LevenbergMarquardtFull;
template<class A> class LevenbergMarquardtReduced;


namespace internal {

template<typename, typename T>
struct has_compute {
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct has_compute<C, Ret(Args...)> {
private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename
        std::is_same<
            decltype(std::declval<T>().compute( std::declval<Args>()... ) ),
            Ret
        >::type;

    template<typename>
    static constexpr std::false_type check(...);
    typedef decltype(check<C>(0)) type;

public:
    static constexpr bool value = type::value;
};


template <template <typename...> class Base, typename Derived>
struct is_derived_from_template {
    using U = typename std::remove_cv<
        typename std::remove_reference<Derived>::type
        >::type;

    template <typename... Args>
        static auto test(Base<Args...>*)
        -> typename std::integral_constant<bool
        , !std::is_same<U, Base<Args...>>::value>;

    static std::false_type test(void*);

    using type = decltype(test(std::declval<U*>()));
};

} // end namespace internal

template <template <typename...> class Base, typename Derived>
using is_derived_from_template
                = typename internal::is_derived_from_template<Base, Derived>::type;

} // end namespace compass

#endif /* FORWARD_DECLARATIONS_HPP */
