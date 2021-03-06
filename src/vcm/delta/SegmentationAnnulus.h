#ifndef SEGMENTATION_ANNULUS_H
#define SEGMENTATION_ANNULUS_H

#include "geometry/DistanceToMeasure.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "geometry/DistanceToMeasureEdge.h"

template<typename ImageFct>
class SegmentationAnnulus {
public:
    typedef typename ImageFct::Value Value;
    typedef typename ImageFct::Point Point;
    typedef typename ImageFct::Domain Domain;
    typedef typename Domain::Space Space;
    typedef typename Space::RealVector RealVector;
    typedef DistanceToMeasureEdge <ImageFct> Distance;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type Container;

public:
    SegmentationAnnulus() = delete;

    SegmentationAnnulus(const Distance &distanceMap,
                        const Container &innerSet);

    SegmentationAnnulus(const SegmentationAnnulus &other);

    ~SegmentationAnnulus();

public:
    Container extractAnnulus();

    Point projection(const Point &p);

private:
    Distance *myDistanceMap;
    Container *myInnerSet;

};


template<typename ImageFct>
SegmentationAnnulus<ImageFct>::
SegmentationAnnulus(const Distance &distanceMap,
                    const Container &innerSet) {
    myDistanceMap = new Distance(distanceMap);
    myInnerSet = new Container(innerSet);
}

template<typename ImageFct>
SegmentationAnnulus<ImageFct>::
SegmentationAnnulus(const SegmentationAnnulus &other) {
    myDistanceMap = new Distance(*other.myDistanceMap);
    myInnerSet = new Container(*other.myInnerSet);
}

template<typename ImageFct>
SegmentationAnnulus<ImageFct>::
~SegmentationAnnulus() {
    if (myDistanceMap) {
        delete myDistanceMap;
        myDistanceMap = 0;
    }
    if (myInnerSet) {
        delete myInnerSet;
        myInnerSet = 0;
    }
}

template<typename ImageFct>
typename SegmentationAnnulus<ImageFct>::Container
SegmentationAnnulus<ImageFct>::
extractAnnulus() {
    Container annulus(myDistanceMap->domain());
    for (const Point &p : *myInnerSet) {
        Point proj = projection(p);
        annulus.insert(proj);
    }
    return annulus;
}

template<typename ImageFct>
typename SegmentationAnnulus<ImageFct>::Point
SegmentationAnnulus<ImageFct>::
projection(const Point &p) {
    RealVector v = myDistanceMap->projection(p);
    Point proj = p + v;
    return proj;
}

#endif
