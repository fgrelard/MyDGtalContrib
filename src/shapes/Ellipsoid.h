#ifndef ELLIPSOID_H
#define ELLIPSOID_H
#include "shapes/Ball.h"
#include "geometry/Distance.h"

template <typename Point>
class Ellipsoid {
public:
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::HyperRectDomain<Space> Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;
    typedef typename Space::RealVector RealVector;
    typedef typename Space::Dimension Dimension;

public:
    Ellipsoid() : myA(0), myB(0), myC(0) {}
    Ellipsoid(const Point& center, double a, double b, double c) : myCenter(center), myA(a), myB(b), myC(c) {}
    Ellipsoid(const Ellipsoid& other) : myCenter(other.myCenter), myA(other.myA), myB(other.myB), myC(other.myC) {}


public:
    inline bool contains(const Point &point) const {
        return std::pow(((point[0] - myCenter[0]) / myA), 2) +std::pow(((point[1] - myCenter[1]) / myB), 2) + std::pow(((point[2] - myCenter[2]) / myC), 2) <= 1.0;
    }

    DigitalSet pointSet() const;

private:
    double myA, myB, myC;
    Point myCenter;
};

template <typename Point>
typename Ellipsoid<Point>::DigitalSet
Ellipsoid<Point>::pointSet() const {
    Ball<Point> ball(myCenter, std::max(myA, std::max(myB, myC)));
    DigitalSet pointBall = ball.pointSet();
    DigitalSet pointEllipsoid(pointBall.domain());
    for (const Point&  p : pointBall) {
        if (contains(p))
            pointEllipsoid.insert(p);
    }
    return pointEllipsoid;
}

#endif
