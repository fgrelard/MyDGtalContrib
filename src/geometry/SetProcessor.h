#ifndef SET_PROCESSOR_H
#define SET_PROCESSOR_H

#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "geometry/Distance.h"
#include "geometry/path/BresenhamAlgorithm.h"
#include "shapes/Ball.h"
#include "geometry/CurveProcessor.h"

template <typename Container>
class SetProcessor {
        BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

public:
        typedef typename Container::value_type Point;
        typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
        typedef typename Space::RealPoint RealPoint;
        typedef typename Space::RealVector RealVector;
        typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
        typedef DGtal::MetricAdjacency<typename Container::Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<typename Container::Space, 3> Adj26;
        typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;

public:
        SetProcessor() = delete;
        SetProcessor(const Container& container)  {
                myContainer = new Container(container);
        }
        SetProcessor(const SetProcessor& other)  {
                myContainer = new Container(*other.myContainer);
        }
        ~SetProcessor() {
                if (myContainer != 0) {
                        delete myContainer;
                        myContainer = 0;
                }
        }

public:
        std::pair<Point, Point> majorAxis();
        double lengthMajorAxis();

        Point closestPointAt(const RealPoint& point);
        std::pair<Point, Point> twoClosestPoints(const Container& other);

        bool sameContainer(const Container& otherContainer);

        Container subSet(const Point& extremity, double radius);
        Container subSet(const Container& constraintNotInSet, double radius);

        std::vector<Container> toConnectedComponents();

        Container intersectionNeighborhoodAt(const Point& p, const Container& other);
        Container intersection(const Container& other);


private:
        Container* myContainer;

};

template <typename Container>
std::pair<typename SetProcessor<Container>::Point, typename SetProcessor<Container>::Point>
SetProcessor<Container>::majorAxis() {
        Point p1, p2;
        L2Metric l2Metric;
        double distanceFarthestPoint = 0;
        for (const Point& p : *myContainer) {
                for (const Point& o : *myContainer) {
                        double currentDistance = l2Metric(p, o);
                        if (l2Metric(p, o) > distanceFarthestPoint) {
                                distanceFarthestPoint = currentDistance;
                                p1 = p;
                                p2 = o;
                        }
                }
        }
        std::pair<Point, Point> pair = std::make_pair(p1, p2);
        return pair;
}


template <typename Container>
double
SetProcessor<Container>::lengthMajorAxis() {
        L2Metric l2Metric;
        std::pair<Point, Point> farthest = majorAxis();
        double distanceFarthestPoint = l2Metric(farthest.first, farthest.second);
        double radius = distanceFarthestPoint / 2.0;
        return radius;
}

template <typename Container>
typename SetProcessor<Container>::Point
SetProcessor<Container>::closestPointAt(const RealPoint& point) {
        if (myContainer->size() == 0) return point;
        return (*std::min_element(myContainer->begin(), myContainer->end(), [&](const Point& p1, const Point& p2) {
                                return (Distance::euclideanDistance((RealPoint)p1, point) < Distance::euclideanDistance((RealPoint)p2, point));
                        }));
}

template <typename Container>
std::pair<typename SetProcessor<Container>::Point, typename SetProcessor<Container>::Point>
SetProcessor<Container>::
twoClosestPoints(const Container& other) {
        L2Metric l2Metric;
        double distance = std::numeric_limits<double>::max();
        std::pair<Point, Point> candidate;
        for (const Point& o : other) {
                Point p = closestPointAt(o);
                double currentDistance = l2Metric(p, o);
                if (currentDistance < distance) {
                        distance = currentDistance;
                        candidate = std::make_pair(p, o);
                }
        }
        return candidate;
}

template <typename Container>
bool
SetProcessor<Container>::sameContainer(const Container & container) {
        if (myContainer->size() != container.size()) return false;
        return (intersection(container).size() == myContainer->size());
}

template <typename Container>
Container
SetProcessor<Container>::
subSet(const Point& extremity, double radius) {
        Ball<Point> ball(extremity, radius);
        Container subSet = ball.intersection(*myContainer);
        if (subSet.size() < 2)
                return *myContainer;
        return subSet;
}

template <typename Container>
Container
SetProcessor<Container>::
subSet(const Container& constraintNotInSet, double radius) {
        Container ep = CurveProcessor<Container>(*myContainer).endPoints();
        Point e;
        for (const Point& p : ep) {
                if (constraintNotInSet.find(p) == constraintNotInSet.end())
                        e = p;
        }
        return subSet(e, radius);
}

template <typename Container>
std::vector<Container> SetProcessor<Container>::toConnectedComponents() {

        Adj26 adj26;
        Adj6 adj6;
        DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
        ObjectType  obj(dt26_6, *myContainer);
        std::vector<ObjectType> cc;
        std::back_insert_iterator<std::vector<ObjectType> > inserter(cc);
        obj.writeComponents(inserter);
        std::vector<Container> ccSet;
        for (const auto& o : cc)
                ccSet.push_back(o.pointSet());
        return ccSet;
}

template <typename Container>
Container
SetProcessor<Container>::
intersectionNeighborhoodAt(const Point& point, const Container& other) {
        Container intersection(myContainer->domain());

        Adj26 adj26;
        Adj6 adj6;
        DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
        ObjectType  obj(dt26_6, *myContainer);

        std::vector<Point> neighbors;
        std::back_insert_iterator< std::vector<Point> > inserter(neighbors);
        obj.writeNeighbors(inserter, point);
        for (const Point& n : neighbors) {
                if (other.find(n) != other.end())
                        intersection.insert(n);
        }
        return intersection;

}

template <typename Container>
Container
SetProcessor<Container>::
intersection(const Container& other) {
        Container intersection(myContainer->domain());
        for (const Point& p : *myContainer) {
                if (other.find(p) != other.end()) {
                        intersection.insert(p);
                }
        }
        return intersection;
}

#endif
