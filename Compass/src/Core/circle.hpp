#ifndef CIRCLE_H
#define CIRCLE_H

#include <eigen-master/Eigen/Dense>

namespace compass {

template<typename T>
class Circle {

public:
    Circle() = default;
    Circle(T a, T b, T radius) : a(a) , b(b), radius(radius) {}
    Circle(const Eigen::Ref<const Eigen::RowVector3<T>>& parameters)
        : a(parameters(0)), b(parameters(1)), radius(parameters(2)) {}

    void setParameters(T a, T b, T radius) {
        this->a = a;
        this->b = b;
        this->radius = radius;
    }

    void setParameters(const Eigen::Ref<const Eigen::RowVector3<T>>& parameters) {
        setParameters(parameters(0), parameters(1), parameters(2));
    }

    void setParameters(T a, T b, T x, T y) {}

    friend std::ostream& operator<<(std::ostream& os, const Circle& circle) {
        return os << "Center: (" << circle.a << "," << circle.b << ")\n" <<
            "Radius: " << circle.radius << std::endl;
    }

    Eigen::Vector3<T> getVector() const {
        return Eigen::Vector3<T>(a, b, radius);
    }

    Eigen::RowVector2<T> getCenter() const {
        return Eigen::RowVector2<T>{a, b};
    }

    inline T getA() const { return a; }
    inline T getB() const { return b; }
    inline T getRadius() const { return radius; }

private:
    T a, b;     // coordinates of the center
    T radius;   // radius of the circle
    T sigma;    // estimate of sigma (standard deviation)
    T grad;     // the norm of the gradient of the objective function
    int i, j;     // the iteration counters (outer and inner, respectively)

};

}
#endif
