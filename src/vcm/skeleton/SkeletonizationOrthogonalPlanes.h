#ifndef SKELETONIZATION_ORTHOGONAL_PLANES_H
#define SKELETONIZATION_ORTHOGONAL_PLANES_H


#include "geometry/DistanceToPointFunctor.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "geometry/MedialAxis.h"
#include "shapes/DigitalPlaneSet.h"
#include "geometry/SetProcessor.h"
#include "geometry/CurveProcessor.h"
#include "geometry/SphericalShellIntersection.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "Statistics.h"

template <typename Container, typename PostProcessing >
class SkeletonizationOrthogonalPlanes {
public:
        typedef typename Container::Space Space;
        typedef typename Space::Point Point;
        typedef typename Space::RealPoint RealPoint;
        typedef typename Space::RealVector RealVector;
        typedef typename DGtal::functors::BallConstantPointFunction<Point, double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Container, KernelFunction> PlaneEstimator;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;
        typedef DigitalPlane<Space> Plane;
        typedef DigitalPlaneSet<Space> PlaneSet;
        typedef typename PlaneEstimator::L2Metric L2Metric;
        typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;

public:
        SkeletonizationOrthogonalPlanes() = delete;
        SkeletonizationOrthogonalPlanes(const Container& setVolume,
                                        bool detectJunctions = true);
        SkeletonizationOrthogonalPlanes(const SkeletonizationOrthogonalPlanes& other);
        ~SkeletonizationOrthogonalPlanes();

public:
        Container skeletonize();

public:
        Point trackNextPoint(const PlaneSet& plane);
        void markPoints(const Point& point);
        void markPoints(const Container& points);
        void markPointsBetweenPlanes(const PlaneSet& currentPlane,
                                     const PlaneSet& previousPlane,
                                     double distanceMax = std::numeric_limits<double>::max());
        Container filterIsolatedPoints(int minSize = 1);
        std::vector<Plane> orientEndPoints();

private:
        Container restrictPlaneSet(const PlaneSet& planeSet, double radius);
        Plane orientPlane(const Plane& undirectedPlane,
                          const Container& endPoints);

        bool isPointThin(const Point& p);
        bool is3ShellPoint(const PlaneSet& p, double radius);

private:
        Container* myVolume;
        PlaneEstimator* myPlaneEstimator;
        bool myJunctionDetection;

//Internals
private:
        Container* mySkeleton;
        Container* myMarkedVertices;
        DTL2* myDT;
        Container* myMedialAxis;
        std::vector<Plane> myPlanes;
};


template <typename Container, typename PostProcessing>
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
SkeletonizationOrthogonalPlanes(const Container& setVolume,
                                bool detectJunctions) {
        L2Metric l2Metric;
        myVolume = new Container(setVolume);

        double r = 5;
        double R = 30;
        KernelFunction chi(1.0, r);
        myPlaneEstimator = new PlaneEstimator(*myVolume, chi, R, r);
        myJunctionDetection = detectJunctions;
        mySkeleton = new Container(myVolume->domain());
        myMarkedVertices = new Container(setVolume.domain());
        myDT = new DTL2(myVolume->domain(), *myVolume, l2Metric);
        MedialAxis<Container> ma(*myVolume);
        myMedialAxis = new Container(ma.pointSet());
}

template <typename Container, typename PostProcessing>
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
SkeletonizationOrthogonalPlanes(const SkeletonizationOrthogonalPlanes& other) {
        myVolume = new Container(*other.myVolume);
        myPlaneEstimator = new PlaneEstimator(*other.myPlaneEstimator);
        myJunctionDetection = other.myJunctionDetection;
        mySkeleton = new Container(*other.mySkeleton);
        myMarkedVertices = new Container(*other.myMarkedVertices);
        myDT = new DTL2(*other.myDT);
        myMedialAxis = new Container(*other.myMedialAxis);
        myPlanes = other.myPlanes;
}

template <typename Container, typename PostProcessing>
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
~SkeletonizationOrthogonalPlanes() {
        if (myVolume) {
                delete myVolume;
                myVolume = 0;
        }
        if (myPlaneEstimator) {
                delete myPlaneEstimator;
                myPlaneEstimator = 0;
        }
        if (mySkeleton) {
                delete mySkeleton;
                mySkeleton = 0;
        }
        if (myMarkedVertices) {
                delete myMarkedVertices;
                myMarkedVertices = 0;
        }
        if (myDT) {
                delete myDT;
                myDT = 0;
        }
        if (myMedialAxis) {
                delete myMedialAxis;
                myMedialAxis = 0;
        }
}

template <typename Container, typename PostProcessing>
Container
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
skeletonize() {
        PlaneSet previous;
        L2Metric l2Metric;
        size_t nbPoints = myVolume->size();
        SetProcessor<Container> procMA(*myMedialAxis);
        Container a3ShellPoints(myVolume->domain());
        Point p = *(max_element(myDT->domain().begin(), myDT->domain().end(), [&](const Point& one,
                                                                const Point& two) {
                                        return ((*myDT)(one) < (*myDT)(two));
                                }));
        double distanceMax = (*myDT)(p) * 2.0;
        DGtal::trace.info() << p << std::endl;
        while (myMarkedVertices->size() < nbPoints) {
                DGtal::trace.progressBar( myMarkedVertices->size(),
                                          nbPoints);
                Point closestPoint = procMA.closestPointAt(p);
                double radius = (*myDT)(closestPoint) + 2.0;
                myPlaneEstimator->setRadius(radius);

                Plane plane = myPlaneEstimator->convergentPlaneAt(p, *myVolume, distanceMax);
                Container planePoints = plane.intersectionWithSetOneCC(*myVolume);
                PlaneSet planeSet(plane, planePoints);
                PlaneSet planeSetG = planeSet;

                markPoints(p);
                markPoints(restrictPlaneSet(planeSet, radius));

                RealPoint centerOfMass = Statistics<Container>(planePoints).extractCenterOfMass();
                if (centerOfMass != RealPoint::zero &&
                    l2Metric(centerOfMass, p) <= sqrt(3)) {
                        if (l2Metric(plane.getCenter(), previous.digitalPlane().getCenter()) <= 2 * sqrt(3)) {
                                //markPointsBetweenPlanes(planeSetG, previous, radius);
                        }
                        Point g = SetProcessor<Container>(planePoints).closestPointAt(centerOfMass);
                        RealVector normal = plane.getPlaneEquation().normal();
                        int connexity = plane.getConnexity();
                        Plane planeG(g, normal, connexity);
                        planeSetG = PlaneSet(planeG, planePoints);

                        myPlanes.push_back(planeG);

                        if (myJunctionDetection && is3ShellPoint(planeSetG, radius))
                                a3ShellPoints.insert(g);

                        if (isPointThin(g))
                                mySkeleton->insert(g);
                        previous = planeSetG;
                }
                p = trackNextPoint(planeSetG);
        }
        Container fillHoles = CurveProcessor<Container>(*mySkeleton).ensureOneCC(*myVolume, sqrt(3), 2 * sqrt(3));
        delete mySkeleton;
        mySkeleton = new Container( fillHoles );


        Container filteredSkeleton = filterIsolatedPoints();
        delete mySkeleton;
        mySkeleton = new Container (filteredSkeleton);

        std::vector<Plane> planesEnd = orientEndPoints();
        myPlanes = planesEnd;
        PostProcessing algo(*mySkeleton, a3ShellPoints, *myVolume, myPlanes);
        Container postProcessedSkeleton = algo.postProcess();

        delete mySkeleton;
        mySkeleton = new Container(postProcessedSkeleton);
        return *mySkeleton;
}

template <typename Container, typename PostProcessing>
typename SkeletonizationOrthogonalPlanes<Container, PostProcessing>::Point
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
trackNextPoint(const PlaneSet& planeSet) {
        Point center =  planeSet.digitalPlane().getCenter();
        RealVector direction = planeSet.digitalPlane().getPlaneEquation().normal();
        Container set = planeSet.pointSet();

        Point current = center;
        double scalar = 1;
        while (current == center ||
               set.find(current) != set.end()) {
                current = center + direction * scalar;
                scalar += 0.5;
        }
        if (myMarkedVertices->find(current) != myMarkedVertices->end()) {
                scalar = 1.0;
                current = center;
                while (current == center ||
                       set.find(current) != set.end()) {
                        current = center - direction * scalar;
                        scalar += 0.5;
                }
                if (myMarkedVertices->find(current) != myMarkedVertices->end()) {
                        current = *(max_element(
                                            myVolume->begin(), myVolume->end(),
                                            [&](const Point& one,
                                                const Point& two) {
                                                    return ((*myDT)(one) < (*myDT)(two) &&
                                                            myMarkedVertices->find(one) == myMarkedVertices->end() &&
                                                            myMarkedVertices->find(two) == myMarkedVertices->end());
                                            }));
                }
        }
        return current;
}

template <typename Container, typename PostProcessing>
void
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
markPoints(const Point& point) {
        myMarkedVertices->insert(point);
}

template <typename Container, typename PostProcessing>
void
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
markPoints(const Container& points) {
        myMarkedVertices->insert(points.begin(), points.end());
}

template <typename Container, typename PostProcessing>
void
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
markPointsBetweenPlanes(const PlaneSet& current,
                        const PlaneSet& previous,
                        double distanceMax) {

        L2Metric l2Metric;
        Container difference(myVolume->domain());
        Plane currentPlane = current.digitalPlane();
        Plane previousPlane = previous.digitalPlane();
        RealVector dirCurrent = currentPlane.getPlaneEquation().normal();
        RealVector dirPrevious = previousPlane.getPlaneEquation().normal();
        Point pCurrent = currentPlane.getCenter();
        Point pPrevious = previousPlane.getCenter();

        if (dirPrevious.dot(dirCurrent) > 0) {
                dirCurrent = -dirCurrent;
        }

        for (const Point& p : *myVolume) {
                if (l2Metric(p, pCurrent) > distanceMax) continue;
                if (currentPlane.isPointAbove(p) &&
                    previousPlane.isPointAbove(p))
                        difference.insert(p);
        }
        myMarkedVertices->insert(difference.begin(), difference.end());
}

template <typename Container, typename PostProcessing>
Container
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
filterIsolatedPoints(int minSize) {
        Container filteredSkeleton(mySkeleton->domain());
        SetProcessor<Container> setProc(*mySkeleton);
        std::vector<Container> components = setProc.toConnectedComponents();
        for (const Container& cc : components) {
                if (cc.size() > minSize) {
                        filteredSkeleton.insert(cc.begin(), cc.end());
                }
        }
        return filteredSkeleton;
}

template <typename Container, typename PostProcessing>
std::vector<typename SkeletonizationOrthogonalPlanes<Container, PostProcessing>::Plane>
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
orientEndPoints() {
        std::vector<Plane> planeEndPoints;
        SetProcessor<Container> setProc(*mySkeleton);
        std::vector<Container> components = setProc.toConnectedComponents();
        L2Metric l2Metric;
        for (const Container& cc : components) {
                if (cc.size() < 2) continue;
                CurveProcessor<Container> curveProc(cc);
                Container localEnd = curveProc.endPoints();
                for (const Point& p : localEnd) {
                        Plane plane = *(std::min_element(
                                                  myPlanes.begin(),
                                                  myPlanes.end(), [&](const Plane& plane, const Plane& otherPlane) {
                                                          return l2Metric(plane.getCenter(), p) < l2Metric(otherPlane.getCenter(), p);
                                                  }));
                        plane = orientPlane(plane, localEnd);
                        planeEndPoints.push_back(plane);
                }
        }
        return planeEndPoints;
}

template <typename Container, typename PostProcessing>
Container
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
restrictPlaneSet(const PlaneSet& planeSet, double radius) {
        L2Metric l2Metric;
        Container set = planeSet.pointSet();
        Point center = planeSet.digitalPlane().getCenter();
        Container toMark(set.domain());
        for (const Point& p : set) {
            if (l2Metric(p, center) <= radius)
                toMark.insert(p);
        }
        return toMark;
}

template <typename Container, typename PostProcessing>
typename SkeletonizationOrthogonalPlanes<Container, PostProcessing>::Plane
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
orientPlane(const Plane& undirectedPlane, const Container& endPoints) {

        L2Metric l2Metric;
        Point center = undirectedPlane.getCenter();
        RealVector normal = undirectedPlane.getPlaneEquation().normal();
        Point candidate = *(std::max_element(endPoints.begin(), endPoints.end(), [&](const Point& one, const Point& two) {
                                return (l2Metric(one, center) < l2Metric(two, center));
                        }));
        RealVector dir = (center - candidate).getNormalized();
        normal = (normal.dot(dir) < 0) ? -normal : normal;
        Plane plane(center, normal, undirectedPlane.getConnexity());
        return plane;
}

template <typename Container, typename PostProcessing>
bool
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
isPointThin(const Point& p) {
        Adj26 adj26;
        Adj6 adj6;
        DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
        ObjectType obj(dt26_6, *mySkeleton);
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
        obj.writeNeighbors(inserter, p);
        if (neighbors.size() < 2) {
                return true;
        }
        return false;

}

template <typename Container, typename PostProcessing>
bool
SkeletonizationOrthogonalPlanes<Container, PostProcessing>::
is3ShellPoint(const PlaneSet& planeSet, double radius) {

        Point center = planeSet.digitalPlane().getCenter();
        RealVector normal = planeSet.digitalPlane().getPlaneEquation().normal();
        int connexity = planeSet.digitalPlane().getConnexity();
        Container set = planeSet.pointSet();
        Container minusSet = Plane(center, -normal, connexity).intersectionWithSetOneCC(*myVolume);
        double radiusCurrent = SetProcessor<Container>(set).lengthMajorAxis() + 2.0;
        double radiusCurrentMinus = SetProcessor<Container>(minusSet).lengthMajorAxis() + 2.0;
        double radiusShell = std::max(4.0, std::max(radiusCurrentMinus, radiusCurrentMinus));
        radiusShell *= 1.2;
        double noise = radiusShell / 2.0;
        SphericalShellIntersection<Container> ssi(*myVolume, center, radiusShell);
        if (ssi.degree(noise) >= 3) {
                return true;
        }
        return false;
}

#endif
