#ifndef COMPASS_STATIC_ASSERT_HPP
#define COMPASS_STATIC_ASSERT_HPP

namespace compass {

    template<bool condition>
    struct static_assertion {};

    template<>
    struct static_assertion<true> {
        enum {
            ALGORITHM_MUST_BE_ALGEBRAIC,
            CLASS_MUST_IMPLEMENT_COMPUTE
        };
    };

    #define COMPASS_STATIC_ASSERT(CONDITION, MSG) \
    if (static_assertion<bool(CONDITION)>::MSG) {}


    #define COMPASS_STATIC_ASSERT_ALGEBRAIC_ONLY(TYPE) \
    COMPASS_STATIC_ASSERT(std::is_base_of<AlgebraicFit<TYPE>, TYPE>::value, \
            ALGORITHM_MUST_BE_ALGEBRAIC)


}

#endif /* COMPASS_STATIC_ASSERT_HPP */
