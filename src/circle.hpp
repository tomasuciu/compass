#pragma once
#ifndef CIRCLE_H
#define CIRCLE_H

#include <iostream>
#include <Eigen/Dense>
#include "util.hpp"

namespace compass {

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>

class Circle {

public:
    Circle() = default;
    Circle(T a, T b, T radius) : a(a), b(b), radius(radius) {}
    Circle(const Eigen::Ref<const Eigen::RowVector3<T>>& parameters) : a(parameters(0)), b(parameters(1)), radius(parameters(2)) {}

    void setParameters(T a, T b, T radius) {
        this->a = a;
        this->b = b;
        this->radius = radius;
    }

    void setParameters(T a, T b, T x, T y) {}

    friend std::ostream& operator<<(std::ostream& os, const Circle& circle) {
        return os << "Center: (" << circle.a << "," << circle.b << ")\n" <<
            "Radius: " << circle.radius << std::endl;

        //TODO: Incorporate additional information into circle printout
//            << "\nSigma: " << circle.sigma << "\nGradient: " << circle.grad << "\nIter: (" << circle.i
//            << "," << circle.j << ")\n";
    }

    // Temporary getter functions
    const Eigen::Vector3<T> getVector() const {
        return Eigen::Vector3<T>(a, b, radius);
    }

    Eigen::RowVector2<T> getCenter() const {
        return Eigen::RowVector2<T>{a, b};
    }

    // TODO implement
    const static Circle<T> centerCircle (const Circle<T>& uncenterd, Eigen::Ref<const Eigen::RowVector2<T>> mean) {
        return Circle<T>(0, 0, 0);
    }

    const inline T getA() const { return a; }
    const inline T getB() const { return b; }
    const inline T getRadius() const { return radius; }

private:
    T a, b;     // coordinates of the center
    T radius;   // radius of the circle
    T sigma;    // estimate of sigma (standard deviation)
    T grad;     // the norm of the gradient of the objective function
    int i, j;     // the iteration counters (outer and inner, respectively)

};

}
#endif
