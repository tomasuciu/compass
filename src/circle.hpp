#pragma once
#ifndef CIRCLE_H
#define CIRCLE_H

#include <iostream>

namespace compass {

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>

class Circle {

public:
    Circle() = default;
    Circle(T a, T b, T radius) : a(a), b(b), radius(radius) {}

    void setParameters(T a, T b, T radius) {
        this->a = a;
        this->b = b;
        this->radius = radius;
    }

    friend std::ostream& operator<<(std::ostream& os, const Circle& circle) {
        return os << "Center: (" << circle.a << "," << circle.b << ")\n" <<
            "Radius: " << circle.radius << std::endl;

        //TODO: Incorporate additional information into circle printout
//            << "\nSigma: " << circle.sigma << "\nGradient: " << circle.grad << "\nIter: (" << circle.i
//            << "," << circle.j << ")\n";
    }

private:
    T a, b;     // coordinates of the center
    T radius;   // radius of the circle
    T sigma;    // estimate of sigma (standard deviation)
    T grad;     // the norm of the gradient of the objective function
    int i, j;     // the iteration counters (outer and inner, respectively)

};

}
#endif
