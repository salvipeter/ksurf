#include <Eigen/LU>

#include "kcurve.hh"

using namespace Transfinite;

// Curvature of a quadratic Bezier curve at its endpoint
static double curvature(const Point3D &a, const Point3D &b, const Point3D &c) {
    auto d1 = (c - b) * 2;
    auto d2 = (a - b * 2 + c) * 2;
    return (d1 ^ d2).norm() / std::pow(d1.norm(), 3);
}
 
KappaCurve::KappaCurve(const Point3D &p1, const Vector3D &n1,
                       const Point3D &p2, const Vector3D &n2) {
    // Compute tangent directions in a common plane
    auto b = n1 ^ n2;
    if (b.norm() < epsilon)
        b = (p2 - p1).normalize() ^ n1;
    auto t1 = n1 ^ b, t2 = n2 ^ b;
    if ((p2 - p1) * t1 < 0)
        t1 *= -1;
    if ((p1 - p2) * t2 < 0)
        t2 *= -1;
    t1.normalize(); t2.normalize();

    PointVector cpts(5);
    cpts[0] = p1;
    cpts[4] = p2;
    auto err = [&](double x) {
        // Compute x1,x2 > 0 such that the vertices of the
        // parabolas are at p1 and p2, respectively
        Eigen::MatrixXd A(5, 5);
        A << -2, 0, t1[0], t1[1], t1[2],
             0, -2, t2[0], t2[1], t2[2],
             t1[0] * (1 - x), t2[0] * x, -1, 0, 0,
             t1[1] * (1 - x), t2[1] * x, 0, -1, 0,
             t1[2] * (1 - x), t2[2] * x, 0, 0, -1;
        Eigen::VectorXd b(5);
        auto q = -(p1 * (1 - x) + p2 * x);
        b << p1 * t1, p2 * t2, q[0], q[1], q[2];
        Eigen::VectorXd sol = A.fullPivLu().solve(b);
        double x1 = sol(0), x2 = sol(1);
        // Fill the remaining control points
        cpts[1] = p1 + t1 * x1;
        cpts[2] = { sol(2), sol(3), sol(4) };
        cpts[3] = p2 + t2 * x2;
        // Compute the G2 error
        auto k1 = curvature(cpts[0], cpts[1], cpts[2]);
        auto k2 = curvature(cpts[4], cpts[3], cpts[2]);
        return k1 - k2;
    }; 

    // Bisection to find best x value in (0,1)
    double xl = epsilon, xh = 1 - epsilon, xm = 0;
    size_t iterations = 20;
    for (size_t i = 0; i < iterations; ++i) {
        auto xm_old = xm;
        xm = (xl + xh) / 2;
        if (xm != 0 && std::abs((xm - xm_old) / xm) < epsilon)
            break;
        auto test = err(xl) * err(xm);
        if (test < 0)
            xh = xm;
        else if (test > 0)
            xl = xm;
        else
            break;
    } 

    // Set up private members
    C1 = BSCurve({cpts[0], cpts[1], cpts[2]});
    C2 = BSCurve({cpts[2], cpts[3], cpts[4]});
    u1 = 0.0; u2 = 1.0; // full curves
    l1 = C1.arcLength(u1, 1);
    l2 = C2.arcLength(0, u2);
    ratio = l1 / (l1 + l2);
}

KappaCurve::KappaCurve(const PointVector &cpts, double u1, double u2) : u1(u1), u2(u2) {
    C1 = BSCurve({cpts[0], cpts[1], cpts[2]});
    C2 = BSCurve({cpts[2], cpts[3], cpts[4]});
    l1 = C1.arcLength(u1, 1);
    l2 = C2.arcLength(0, u2);
    ratio = l1 / (l1 + l2);
}

KappaCurve::~KappaCurve() {
}

std::pair<bool, double> KappaCurve::parameter(double u) const {
    if (u < ratio)
        return { true, u / ratio * (1 - u1) + u1 };
    return { false, (u - ratio) / (1 - ratio) * u2 };
}

Point3D KappaCurve::eval(double u) const {
    auto [first, t] = parameter(u);
    if (first)
        return C1.eval(t);
    return C2.eval(t);
}

Point3D KappaCurve::eval(double u, size_t nr_der, VectorVector &der) const {
    auto [first, t] = parameter(u);
    if (first)
        return C1.eval(t, nr_der, der);
    return C2.eval(t, nr_der, der);
}

void KappaCurve::reverse() {
    C1.reverse(); C1.normalize();
    C2.reverse(); C2.normalize();
    std::swap(C1, C2);
    std::swap(u1, u2);
    u1 = 1 - u1;
    u2 = 1 - u2;
    std::swap(l1, l2);
    ratio = 1 - ratio;
}

double KappaCurve::arcLength(double from, double to) const {
    if (from >= to)
        return 0.0;
    auto [first1, t1] = parameter(from);
    auto [first2, t2] = parameter(to);
    if (first2)
        return C1.arcLength(t1, t2);
    if (!first1)
        return C2.arcLength(t1, t2);
    return C1.arcLength(t1, 1) + C2.arcLength(0, t2);
}
