#ifndef BALL_H
#define BALL_H

#include <vector>
#include <set>
#include <iostream>
#include "DGtal/base/Common.h"
#include "geometry/Distance.h"
#include "geometry/PointUtil.h"

template <typename Point>
class Ball {
public:
	typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
	typedef DGtal::HyperRectDomain< Space > Domain;
    typedef typename DGtal::DigitalSetSelector< Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS >::Type DigitalSet;
	typedef typename Space::RealVector RealVector;
public:
	Ball() : myRadius(0.0), myCenter({0,0,0}) {}
	Ball(const Point& center, double radius) : myCenter(center), myRadius(radius) {}
	Ball(const Ball& other) : myCenter(other.myCenter), myRadius(other.myRadius) {}

public:
	bool contains(const Point& point) const {return Distance::euclideanDistance(point, myCenter) <= myRadius;}
	DigitalSet intersection(const DigitalSet& setPoint);
	DigitalSet surfaceIntersection(const DigitalSet& setSurface);
    DigitalSet pointSet() const;
	DigitalSet pointsInHalfBall(const RealVector& normal = RealVector(0,1,0)) const;
	DigitalSet pointsSurfaceBall() const;

public:
	Point getCenter() const {return myCenter;}
	double getRadius() const {return myRadius;}
public:
	bool operator!=(const Ball & other) const {return (myCenter != other.myCenter || myRadius != other.myRadius);}
private:
	Point myCenter;
	double myRadius;
};

template <typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::intersection(const DigitalSet& setPoint) {
	DigitalSet intersection(setPoint.domain());
	for (auto it = setPoint.begin(), ite = setPoint.end(); it != ite; ++it) {
		double distance = Distance::euclideanDistance(*it, myCenter);
		if (distance <= myRadius) {
			intersection.insert(*it);
		}
	}
	return intersection;
}

template <typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::surfaceIntersection(const DigitalSet& setSurface) {
    DigitalSet intersection(setSurface.domain());
	for (auto it = setSurface.begin(), ite = setSurface.end(); it != ite; ++it) {
		double distance = Distance::euclideanDistance(*it, myCenter);
		if (distance >= myRadius-1 && distance <= myRadius) {
			intersection.insert(*it);
		}
	}
	return intersection;
}

template <typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::pointSet() const {
    Point lower(-myRadius + myCenter[0], -myRadius + myCenter[1], -myRadius + myCenter[2]);
	Point upper(myRadius + myCenter[0] + 1, myRadius + myCenter[1] + 1, myRadius + myCenter[2] + 1);
    Domain domain(lower, upper);
	DigitalSet points(domain);
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				if (contains(p)) {
					points.insert(p);
				}
			}
		}
	}
	return points;
}



template <typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::pointsInHalfBall(const RealVector& normal) const {
	std::set<Point> points;
	double d = myCenter[0] * normal[0] + myCenter[1] * normal[1] + myCenter[2] * normal[2];
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				double eq = p[0] * normal[0] + p[1] * normal[1] + p[2] * normal[2] - d;
				if (contains(p) && eq > 0) {
					points.insert(p);
				}
			}
		}
	}
	Domain domain = PointUtil::computeBoundingBox<Domain>(points);
    DigitalSet aSet(domain);
	aSet.insert(points.begin(), points.end());
	return aSet;
}

template <typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::pointsSurfaceBall() const {
	std::set<Point> points;
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				double distance = Distance::euclideanDistance(p, myCenter);
				if (distance >= myRadius-1 && distance <= myRadius) {
					points.insert(p);
				}
			}
		}
	}
	Domain domain = PointUtil::computeBoundingBox<Domain>(points);
    DigitalSet aSet(domain);
	aSet.insert(points.begin(), points.end());
	return aSet;
}
#endif
