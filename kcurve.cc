#include "kcurve.hh"

using namespace Transfinite;
 
KappaCurve::KappaCurve(const Point3D &p1, const Vector3D &n1,
                       const Point3D &p2, const Vector3D &n2) {
    // Compute tangent directions
    auto d1 = p1 + (p2 - p1) / 3, d2 = p2 + (p1 - p2) / 3;
    d1 -= n1 * ((d1 - p1) * n1); d2 -= n2 * ((d2 - p2) * n2);
    auto t1 = (d1 - p1).normalize(), t2 = (d2 - p2).normalize();

    // Set cpts, u1, u2  based on  p1, p2, t1, t2
    PointVector cpts;
    u1 = 0.0; u2 = 1.0; // full curves

    // Set up private members
    C1 = BSCurve({cpts[0], cpts[1], cpts[2]});
    C2 = BSCurve({cpts[2], cpts[3], cpts[4]});
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
