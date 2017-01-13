#ifndef JUNCTION_PROCESSING_SKELETON_H
#define JUNCTION_PROCESSING_SKELETON_H

#include <vector>
#include "vcm/skeleton/PostProcessingSkeleton.h"
#include "shapes/DigitalPlane.h"
#include "geometry/TruePredicate.h"
#include "geometry/SetProcessor.h"


template <typename Container, typename Predicate = TruePredicate<typename Container::Space> >
class JunctionProcessingSkeleton {

public:
    typedef typename Container::Space Space;
    typedef typename Container::Point Point;
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef typename Space::RealVector RealVector;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;
    typedef DigitalPlane<Space> Plane;

public:
    JunctionProcessingSkeleton(const Container& skeletonPoints,
                               const Container& a3ShellPoints,
                               const Container& setVolume,
                               const std::vector<Plane>& planesEndPoints);
    ~JunctionProcessingSkeleton();

public:
    Container postProcess();

private:
    Container shellPointsToJunctionAreas();
    Container dilate(const Container& toDilate);
    std::vector<Container> groupsOfPointSameJunction();
    Point referencePointWithDifferenceNormal(const Container& points);

private:
    Container* mySkeleton;
    Container* my3ShellPoints;
    Container* myVolume;
    Container* myJunctionArea;
    std::vector<Plane>* myPlanes;
};


template <typename Container, typename Predicate>
JunctionProcessingSkeleton<Container, Predicate>::
JunctionProcessingSkeleton(const Container& skeletonPoints,
                           const Container& a3ShellPoints,
                           const Container& setVolume,
                           const std::vector<Plane>& planesEndPoints) {
        mySkeleton = new Container(skeletonPoints);
        my3ShellPoints = new Container(a3ShellPoints);
        myVolume = new Container(setVolume);
        Container firstJunctionArea = new Container( shellPointsToJunctionAreas() );

        delete my3ShellPoints;
        my3ShellPoints = 0;

        my3ShellPoints = new Container( dilate(firstJunctionArea) );
        myJunctionArea = new Container( shellPointsToJunctionAreas() );
        myPlanes = new std::vector<Plane>(planesEndPoints);
}

template <typename Container, typename Predicate>
JunctionProcessingSkeleton<Container, Predicate>::
~JunctionProcessingSkeleton() {
        if (mySkeleton) {
                delete mySkeleton;
                mySkeleton = 0;
        }
        if (my3ShellPoints) {
                delete my3ShellPoints;
                my3ShellPoints = 0;
        }
        if (myVolume) {
                delete myVolume;
                myVolume = 0;
        }
        if (myJunctionArea) {
                delete myJunctionArea;
                myJunctionArea = 0;
        }
        if (myPlanes) {
                delete myPlanes;
                myPlanes = 0;
        }
}



template <typename Container, typename Predicate>
Container
JunctionProcessingSkeleton<Container, Predicate>::
shellPointsToJunctionAreas() {
        L2Metric l2Metric;
        Container shellAreas(myVolume.domain());
        for (const Point& p : *myVolume) {
                Point closestPoint = *min_element(mySkeleton->begin(), mySkeleton->end(), [&](const Point& one, const Point& two) {
                                return l2Metric(one, p) < l2Metric(two, p);
                        });
                if (my3ShellPoints->find(closestPoint) != my3ShellPoints->end())
                        shellAreas.insert(p);
        }
        return shellAreas;
}


template <typename Container, typename Predicate>
Container
JunctionProcessingSkeleton<Container, Predicate>::
dilate(const Container& toDilate) {
    Container dilatedJunction(toDilate.domain());

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
    ObjectType obj(dt26_6, *myVolume);

    for (const Point& s : *mySkeleton) {
            std::vector<Point> neighbors;
            std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
            obj.writeNeighbors(inserter, s);
            for (const Point& n : neighbors) {
                    if (toDilate.find(n) != toDilate.end())
                             toDilate.insert(s);
            }

    }
    return toDilate;
}


template <typename Container, typename Predicate>
std::vector<Container>
JunctionProcessingSkeleton<Container, Predicate>::
groupsOfPointSameJunction() {

    std::vector< Container > groups;
    Predicate predicate(myPlanes);

    std::vector<Container> ccJunction = SetProcessor<Container>(myJunctionArea).toConnectedComponents();

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
    ObjectType obj(dt26_6, *mySkeleton);

    std::map<Point, int> labelMap;
    for (const Plane& plane : *myPlanes) {
        Point currentPoint = plane.getCenter();
        double currentDistance = std::numeric_limits<double>::max();
        int indexJunction;
        int i = -1;
        for (const Container& junction : ccJunction) {
            i++;
            Point closest = SetProcessor<Container>(junction).closestPointAt(currentPoint);

            bool add = false;
            for (const Point& p : junction) {
                add |= predicate(p);
            }

            double distance = l2Metric(currentPoint, closest);
            if (currentDistance > distance && add) {
                currentDistance = distance;
                indexJunction = i;
            }
        }
        if (indexJunction >= 0)
            labelMap[currentPoint] = i;
    }

    Container emptyContainer(myVolume->domain());
    groups.resize(ccJunction.size(), emptyContainer);
    for (const std::pair<Point, int>& pointToJunctionPos : labelMap) {
        groups[pointToJunctionPos.second].insert(pointToJunctionPos.first);
    }
    return groups;
}

template <typename Container, typename Predicate>
typename JunctionProcessingSkeleton<Container, Predicate>::Point
JunctionProcessingSkeleton<Container, Predicate>::
referencePointWithDifferenceNormal(const Container& points) {
    int referenceNumber = std::numeric_limits<int>::min();
    Point cand = *(points.begin());
    for (const Point& p : points) {
        auto itPlane = find_if(myPlanes->begin(), myPlanes.end(),
                               [&](const Plane& plane) {
                                   return (plane.getCenter() == p);
                               });
        if (itPlane == myPlanes->end()) continue;
        RealVector n = itPlane->getPlaneEquation().normal();
        int cpt = 0;
        for (const Point& p2 : points) {
            if (p == p2) continue;
            auto itPlane2 = find_if(myPlanes->begin(), myPlanes.end(),
                                   [&](const Plane& plane) {
                                       return (plane.getCenter() == p2);
                                   });
            if (itPlane2 == myPlanes->end()) continue;
            RealVector n2 = itPlane->getPlaneEquation().normal();
            if (n.dot(n2) < 0)
                cpt++;
        }
        if (cpt > referenceNumber) {
            referenceNumber = cpt;
            cand = p;
        }
    }
    return cand;
}

#endif
