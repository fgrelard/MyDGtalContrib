#ifndef SPHERICAL_SHELL_INTERSECTION_H
#define SPHERICAL_SHELL_INTERSECTION_H

#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/src/DGtal/topology/Object.h"

template <typename Container>
class SphericalShellIntersection {
        BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

public:
        typedef typename Container::Point Point;
        typedef typename Container::Space Space;
        typedef typename Space::RealVector RealVector;
        typedef typename DGtal::HyperRectDomain< Space > Domain;
        typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;
public:
        SphericalShellIntersection() = delete;
        SphericalShellIntersection(const Point& center, const Container& setVolume, double radiusInnerBall, double radiusOuterBall, const RealVector& dirVector = RealVector::zero);
        SphericalShellIntersection(const SphericalShellIntersection& other);
        ~SphericalShellIntersection();

private:
        Container* myContainer;
        Point myCenter;
        double myRadiusInner;
        double myRadiusOuter;
        RealVector myDirectionVector;
}

#endif
