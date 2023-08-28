#pragma once

#include "curves.hh"

class KappaCurve : public Transfinite::Curve {
public:
    KappaCurve(const Transfinite::Point3D &p1,
               const Transfinite::Vector3D &n1,
               const Transfinite::Point3D &p2,
               const Transfinite::Vector3D &n2);
    // cpts contains 5 control points, defining 2 quadratic Bezier curves:
    // C1 : cpts[0], cpts[1], cpts[2]
    // C2 : cpts[2], cpts[3], cpts[4]
    // The kappa-curve is the conjunction of C1 and C2
    // in the respective intervals [u1...1] and [0...u2]
    KappaCurve(const Transfinite::PointVector &cpts, double u1, double u2);
    ~KappaCurve();
    Transfinite::Point3D eval(double u) const override;
    Transfinite::Point3D eval(double u, size_t nr_der,
            Transfinite::VectorVector &der) const override;
    void reverse() override;
    double arcLength(double from, double to) const override;
private:
    std::pair<bool, double> parameter(double u) const;

    Transfinite::BSCurve C1, C2;
    double u1, u2, l1, l2, ratio;
};


