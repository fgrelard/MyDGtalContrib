#ifndef DISTANCE_TO_MEASURE_EDGE_H
#define DISTANCE_TO_MEASURE_EDGE_H

#include "DGtal/base/Trace.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "geometry/DistanceToPointFunctor.h"
#include "PointUtil.h"
#include "CurveProcessor.h"


template<typename ImageFct>
class DistanceToMeasureEdge {
public:
    typedef typename ImageFct::Value Value;
    typedef typename ImageFct::Point Point;
    typedef typename ImageFct::Domain Domain;
    typedef typename Domain::Space Space;
    typedef typename Space::Dimension Dimension;
    typedef typename Space::RealVector RealVector;
    typedef typename Point::Scalar Scalar;
    typedef DGtal::MetricAdjacency<Space, 1> Adjacency;

public:

    DistanceToMeasureEdge(Value m0, const ImageFct &measure, Value rmax = 10.0, bool initialize = true)
            : myMass(m0), myMeasure(measure), myDistance2(myMeasure.domain()),
              myR2Max(rmax * rmax) {
        if (initialize)
            init();
    }


    DistanceToMeasureEdge(const DistanceToMeasureEdge &other) : myMass(other.myMass), myMeasure(other.myMeasure),
                                                                myDistance2(other.myDistance2), myR2Max(other.myR2Max) {
    }

    virtual void init() {
        double nb = myDistance2.domain().size();
        unsigned int i = 0;

        for (typename Domain::ConstIterator it = myDistance2.domain().begin(),
                     itE = myDistance2.domain().end(); it != itE; ++it, ++i) {
            DGtal::trace.progressBar(i, nb);
            myDistance2.setValue(*it, computeDistance2(*it));
        }
    }

    inline const ImageFct &measure() const {
        return myMeasure;
    }

    inline Value operator()(const Point &p) const {
        return distance(p);
    }

    Value distance(const Point &p) const {
        return sqrt(distance2(p));
    }


    Value distance2(const Point &p) const {
        return myDistance2(p);
    }

    Value safeDistance2(const Point &p) const {
        if (myDistance2.domain().isInside(p))
            return myDistance2(p);
        else return myDistance2(box(p));
    }

    Point box(const Point &p) const {
        Point q = p.sup(myDistance2.domain().lowerBound());
        return q.inf(myDistance2.domain().upperBound());
    }

    virtual RealVector projection(const Point &p) const {
        std::vector<Point> neighborsP;
        std::back_insert_iterator<std::vector<Point> > outIterator(neighborsP);
        Adjacency::writeNeighbors(outIterator, p);

        typedef typename std::vector<Point>::iterator Iterator;
        Value distance_center = distance2(p);
        RealVector vectorToReturn;
        for (Iterator it = neighborsP.begin(), ite = neighborsP.end();
             it != ite; ++it) {
            Point n = *it;
            Value distance = (myDistance2.domain().isInside(n)) ? distance2(n) : distance_center;
            Point diff = p - n;
            Point otherPoint = p + diff;
            Value otherDistance = (myDistance2.domain().isInside(otherPoint)) ? distance2(otherPoint)
                                                                              : distance_center;
            auto valMax = std::max_element(diff.begin(), diff.end(), [&](const Scalar &one, const Scalar &two) {
                return std::abs(one) < std::abs(two);
            });
            int d = valMax - diff.begin();

            if (otherPoint[d] > n[d]) {
                Value tmpDistance = otherDistance;
                otherDistance = distance;
                distance = tmpDistance;
            }

            //Necessary to avoid incorrect vector direction for pixels of max distance near pixels of min distance (background)
            if (std::abs(std::sqrt(distance) - std::sqrt(distance_center)) > Point::dimension) {
                distance = distance_center;
            }
            if (std::abs(std::sqrt(otherDistance) - std::sqrt(distance_center)) > Point::dimension) {
                otherDistance = distance_center;
            }

            vectorToReturn[d] = (std::abs(distance - distance_center) >=
                                 std::abs(distance_center - otherDistance)) ? -(distance - distance_center) /
                                                                              2.0 :
                                -(distance_center - otherDistance) / 2.0;


        }
        return vectorToReturn;
    }

    virtual RealVector projectionDistance(const Point &p) const {
        RealVector proj = projection(p);
        return proj.getNormalized() * sqrt(myDistance2(p));
    }

    virtual Value computeDistance2(const Point &p) {
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> Distance;
        typedef DistanceToPointFunctor<Distance> DistanceToPoint;

        typedef DGtal::DistanceBreadthFirstVisitor<Adjacency, DistanceToPoint, std::set<Point> >
                DistanceVisitor;
        typedef typename DistanceVisitor::Node MyNode;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology<Adj26, Adj6> Topology;
        typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;
        typedef DGtal::Object<Topology, DigitalSet> ObjectType;

        Adj26 adj26;
        Adj6 adj6;
        Topology dt26_6(adj26, adj6, DGtal::JORDAN_DT);

        Value m = DGtal::NumberTraits<Value>::ZERO;
        Value d2 = DGtal::NumberTraits<Value>::ZERO;
        Adjacency graph;
        Distance l2;
        DistanceToPoint d2pfct(l2, p);
        DistanceVisitor visitor(graph, d2pfct, p);

        Value last = d2pfct(p);
        MyNode node;
        Value firstMass = myMeasure(p);
        DGtal::Statistic<Value> stat(true);
        stat.addValue(firstMass);
        while (!visitor.finished()) {
            node = visitor.current();
            std::vector<MyNode> vec;
            visitor.getCurrentLayer(vec);

            for (const MyNode &n : vec) {
                if (!myMeasure.domain().isInside(n.first)) continue;
                double currentColor = myMeasure(n.first);
                stat.addValue(currentColor);
            }
            firstMass = stat.median();
            m = DGtal::NumberTraits<Value>::ZERO;
            for (const Value &v : stat) {
                m += v - firstMass;
            }

            if (node.second >= std::sqrt(myR2Max) / 2 && m < 0) {
                node.second = myR2Max;
                break;
            }
            if (m >= myMass) {
                break;
            }
            if (node.second > std::sqrt(myR2Max)) {
                d2 = m * myR2Max;
                break;
            }
            visitor.expandLayer();
        }
        if (m == DGtal::NumberTraits<Value>::ZERO)
            return 0.0;
        if (node.second == myR2Max)
            return 0.0;
        return node.second * node.second;
    }

    Domain domain() const { return myMeasure.domain(); }

protected:
    Value myMass;
    const ImageFct &myMeasure;
    ImageFct myDistance2;
    Value myR2Max;
};

#endif
