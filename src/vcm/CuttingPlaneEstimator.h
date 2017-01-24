#ifndef CUTTING_PLANE_ESTIMATOR_H
#define CUTTING_PLANE_ESTIMATOR_H

#include <vector>
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "geometry/CurveProcessor.h"

template <typename PlaneEstimator>
class CuttingPlaneEstimator {
public:
        typedef typename PlaneEstimator::Plane Plane;
        typedef typename PlaneEstimator::Container Container;
        typedef typename PlaneEstimator::Point Point;
        typedef typename PlaneEstimator::RealVector RealVector;
        typedef typename PlaneEstimator::Domain Domain;
        typedef typename PlaneEstimator::Space Space;
        typedef typename PlaneEstimator::L2Metric L2Metric;
        typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;

public:
        CuttingPlaneEstimator();
        CuttingPlaneEstimator(const std::vector<Container>& branches,
                      const PlaneEstimator& planeEstimator,
                      const Container& setVolume,
                      const DTL2& aDT) ;
        CuttingPlaneEstimator(const CuttingPlaneEstimator& other);
        ~CuttingPlaneEstimator();

public:
        std::vector<Plane> cuttingPlanes();

public:
        std::vector<Plane> computePlanes(const std::vector<Point>& orientedEdge);
        Plane cuttingPlane(const std::vector<Plane>& planes);
        Plane referencePlane(const Point& referencePoint);
        std::vector<Plane> filteredPlanes(const std::vector<Plane>& cuttingPlanes,
                                          const std::vector<Plane>& endPlanes,
                                          const Plane& referencePlane);
        Point extractBranchingPoint();

private:
        std::vector<Container>* myBranches;
        PlaneEstimator* myPlaneEstimator;
        Container* myVolume;
        DTL2* myDT;

};

template <typename PlaneEstimator>
CuttingPlaneEstimator<PlaneEstimator>::
CuttingPlaneEstimator() {
        myBranches = 0;
        myPlaneEstimator = 0;
        myVolume = 0;
        myDT = 0;
}

template <typename PlaneEstimator>
CuttingPlaneEstimator<PlaneEstimator>::
CuttingPlaneEstimator(const std::vector<Container>& branches,
              const PlaneEstimator& planeEstimator,
              const Container& setVolume,
              const DTL2& aDT) {
        myBranches = new std::vector<Container>( branches );
        myPlaneEstimator = new PlaneEstimator( planeEstimator );
        myVolume = new Container( setVolume );
        myDT = new DTL2( aDT );
}

template <typename PlaneEstimator>
CuttingPlaneEstimator<PlaneEstimator>::
CuttingPlaneEstimator( const CuttingPlaneEstimator& other) {
        myBranches = new std::vector<Container>( *other.myBranches );
        myPlaneEstimator = new PlaneEstimator( *other.myPlaneEstimator );
        myVolume = new Container( *other.myContainer );
        myDT = new DTL2 ( *other.myDT );
}

template <typename PlaneEstimator>
CuttingPlaneEstimator<PlaneEstimator>::
~CuttingPlaneEstimator() {
        if (myBranches) {
                delete myBranches;
                myBranches = 0;
        }
        if (myPlaneEstimator) {
                delete myPlaneEstimator;
                myPlaneEstimator = 0;
        }
        if (myVolume) {
                delete myVolume;
                myVolume = 0;
        }
        if (myDT) {
                delete myDT;
                myDT = 0;
        }
}

template <typename PlaneEstimator>
std::vector<typename CuttingPlaneEstimator<PlaneEstimator>::Plane>
CuttingPlaneEstimator<PlaneEstimator>::
cuttingPlanes() {
        std::vector<Plane> cuttingPlanes;
        std::vector<Plane> endPlanes;
        Point b = extractBranchingPoint();

        for (const Container& branch : *myBranches) {
                std::vector<Point> orderedBranch = CurveProcessor<Container>(branch).convertToOrderedCurve(b);
                std::vector<Plane> planes = computePlanes(orderedBranch);
                Plane cutting = cuttingPlane(planes);
                Plane end = *(planes.rbegin());

                cuttingPlanes.push_back(cutting);
                endPlanes.push_back(end);
        }
        Plane planeReference = referencePlane(b);
        cuttingPlanes = filteredPlanes(cuttingPlanes, endPlanes, planeReference);
        return cuttingPlanes;

}

template <typename PlaneEstimator>
std::vector<typename CuttingPlaneEstimator<PlaneEstimator>::Plane>
CuttingPlaneEstimator<PlaneEstimator>::
computePlanes(const std::vector<Point>& orientedEdge) {
        std::vector<Plane> planes;
        for (const Point& p : orientedEdge) {
                double radius = (*myDT)(p) + 2;
                myPlaneEstimator->setRadius(radius);
                Plane plane = myPlaneEstimator->convergentPlaneAt(p, *myVolume, radius*2);
                planes.push_back(plane);

        }
        return planes;
}

template <typename PlaneEstimator>
typename CuttingPlaneEstimator<PlaneEstimator>::Plane
CuttingPlaneEstimator<PlaneEstimator>::
cuttingPlane(const std::vector<Plane>& planes) {
        typedef DGtal::MetricAdjacency<Space, 3> MAdj;

        Plane candidate;
        double previousFactor = 0;
        for (const Plane& plane : planes) {
                Point p = plane.getCenter();
                Container current = plane.intersectionWithSetOneCC(*myVolume);
                double currentValue = current.size();
                std::vector<Point> neighbors;
                std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
                MAdj::writeNeighbors(inserter, p);
                for (const Point& n : neighbors) {
                        auto nInMap = find_if(planes.begin(), planes.end(), [&](const Plane& planeN) {
                                        return planeN.getCenter() == n;
                                });
                        if (nInMap != planes.end()) {
                                double valueNeighbor = nInMap->intersectionWithSetOneCC(*myVolume).size();
                                double factor = valueNeighbor / currentValue;
                                if (factor > previousFactor ) {
                                        previousFactor = factor;
                                        candidate = plane;
                                }
                        }
                }
        }
        return candidate;

}

template <typename PlaneEstimator>
typename CuttingPlaneEstimator<PlaneEstimator>::Plane
CuttingPlaneEstimator<PlaneEstimator>::
referencePlane(const Point& referencePoint) {
        double radius = (*myDT)(referencePoint) + 2;
        myPlaneEstimator->setRadius(radius);
        return myPlaneEstimator->convergentPlaneAt(referencePoint, *myVolume, radius * 5);
}

template <typename PlaneEstimator>
std::vector<typename CuttingPlaneEstimator<PlaneEstimator>::Plane>
CuttingPlaneEstimator<PlaneEstimator>::
filteredPlanes(const std::vector<Plane>& cuttingPlanes, const std::vector<Plane>& endPlanes, const Plane& referencePlane) {

        std::vector<Plane> first, second;

        for (int i = 0, e  = endPlanes.size(); i < e; i++) {
                Plane endPlane = endPlanes[i];
                Plane initPlane = cuttingPlanes[i];
                if (referencePlane.contains(endPlane.getCenter())) continue;
                if (referencePlane.isPointAbove(endPlane.getCenter())) {
                        first.push_back(initPlane);
                }
                else {
                        second.push_back(initPlane);
                }
        }
        return (second.size() > first.size()) ? second : first;
}

template <typename PlaneEstimator>
typename CuttingPlaneEstimator<PlaneEstimator>::Point
CuttingPlaneEstimator<PlaneEstimator>::
extractBranchingPoint() {
        for (int i = 0; i < myBranches->size(); i++) {
                Container setCurrent = (*myBranches)[i];
                for (int j = i + 1; j < myBranches->size(); j++) {
                        Container setOther = (*myBranches)[j];
                        Container intersection = SetProcessor<Container>(setCurrent).intersection(setOther);
                        if (intersection.size() == 1)
                                return *intersection.begin();
                }
        }
        return Point::zero;
}


#endif
