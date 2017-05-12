#ifndef DISTANCE_TO_MEASURE_EDGE_LOCAL_H
#define DISTANCE_TO_MEASURE_EDGE_LOCAL_H

#include "DGtal/base/Trace.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "geometry/DistanceToPointFunctor.h"
#include "DGtal/math/Statistic.h"
#include "PointUtil.h"
#include "CurveProcessor.h"
#include <omp.h>

template<typename ImageFct>
class DistanceToMeasureEdgeLocal {
public:
    typedef typename ImageFct::Value Value;
    typedef typename ImageFct::Point Point;
    typedef typename ImageFct::Domain Domain;
    typedef typename Domain::Space Space;
    typedef typename Space::Dimension Dimension;
    typedef typename Space::RealVector RealVector;
    typedef typename Point::Scalar Scalar;
    typedef DGtal::MetricAdjacency<Space, 1> Adjacency;
    typedef DGtal::MetricAdjacency<Space, 3> Adjacency26;

public:
    DistanceToMeasureEdgeLocal() : myMeasure(Domain(Point::zero, Point::zero)), myDistance2(Domain(Point::zero, Point::zero)) {}

    DistanceToMeasureEdgeLocal(const ImageFct& distance2) : myDistance2(distance2), myMeasure(distance2.domain()) {}

    DistanceToMeasureEdgeLocal(Value m0, const ImageFct &measure, Value rmax = 10.0, Value mask = 1.1, bool initialize =
    true)
            : myMass(m0), myMeasure(measure), myDistance2(myMeasure.domain()),
              myR2Max(rmax * rmax), myMask(mask) {
        if (initialize)
            init();
    }


    DistanceToMeasureEdgeLocal(const DistanceToMeasureEdgeLocal &other) : myMass(other.myMass), myMeasure(other.myMeasure),
                                                                myDistance2(other.myDistance2), myR2Max(other.myR2Max),
                                                                myMask(other.myMask) {
    }

    virtual void init() {
        int nb = myDistance2.domain().size();
        Domain domain = myMeasure.domain();
        int i = 0;
        for (typename Domain::Iterator it = domain.begin(), ite = domain.end(); it != ite; ++it, ++i) {
            DGtal::trace.progressBar(i, nb);
            Point current = *it;
            bool isMasked = (myMeasure(current) == myMask);

            if (!isMasked)
                myDistance2.setValue(current, computeDistance2(current));
            else
                myDistance2.setValue(current, 0);
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
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, 1> Distance;
        typedef DistanceToPointFunctor<Distance> DistanceToPoint;

        typedef DGtal::DistanceBreadthFirstVisitor<Adjacency, DistanceToPoint, std::set<Point> >
                DistanceVisitor;
        typedef typename DistanceVisitor::Node MyNode;
        typedef DGtal::DigitalSetBySTLVector<Domain> Set;

        Value m = DGtal::NumberTraits<Value>::ZERO;
        Adjacency graph;
        Distance l2;
        DistanceToPoint d2pfct(l2, p);
        DistanceVisitor visitor(graph, d2pfct, p);
        Value last = d2pfct(p);
        MyNode node;
        Value firstMass = myMeasure(p);
        DGtal::Statistic<Value> stat(true), statNext(true);
        stat.addValue(firstMass);

        Value previousMean = DGtal::NumberTraits<Value>::ZERO;
        Value currentMean = firstMass;
        while (!visitor.finished()) {
            node = visitor.current();

            std::vector<MyNode> vec;
            visitor.getCurrentLayer(vec);

            std::vector<Point> neighbors;
            std::back_insert_iterator<std::vector<Point> > inserter(neighbors);

            for (const MyNode &n : vec) {
                if (!myMeasure.domain().isInside(n.first)) continue;
                double currentColor = myMeasure(n.first);
                stat.addValue(currentColor);
            }

            currentMean = stat.mean() + std::sqrt(stat.variance());
            for (const MyNode &n : vec) {
                if (!myMeasure.domain().isInside(n.first)) continue;
                double currentColor = myMeasure(n.first);
                if (currentColor > currentMean) {
                    Adjacency26::writeNeighbors(inserter, n.first);
                }
            }


            visitor.expandLayer();

            std::vector<MyNode> vecNext;
            visitor.getCurrentLayer(vecNext);


            Set aSet(myMeasure.domain());
            for (const MyNode &n : vecNext) {
                if (!myMeasure.domain().isInside(n.first)) continue;
                if (std::find(neighbors.begin(), neighbors.end(), n.first) == neighbors.end()) continue;
                statNext = stat;
                statNext.addValue(myMeasure(n.first));
                if (statNext.variance() > stat.variance() && statNext.mean() > stat.mean())
                    aSet.insert(n.first);
            }

            SetProcessor<Set> processor(aSet);
            std::vector<Set> cc = processor.toConnectedComponents();
            for (const Set& s : cc) {
                if (s.size() >= 2)
                    m++;
            }

            if (m > myMass)
                break;
            if (node.second * node.second > myR2Max)
                break;

        }
        // if (m == DGtal::NumberTraits<Value>::ZERO)
        //     return DGtal::NumberTraits<Value>::ZERO;
        // if (currentMean < previousMean)
        //     return DGtal::NumberTraits<Value>::ZERO;

        return node.second * node.second;
    }

    Domain domain() const { return myMeasure.domain(); }

    DistanceToMeasureEdgeLocal& operator=(const DistanceToMeasureEdgeLocal& other) {
        myMass = other.myMass;
        myMask = other.myMask;
        myMeasure = other.myMeasure;
        myDistance2 = other.myDistance2;
        myR2Max = other.myR2Max;
        return *this;
    }

protected:
    Value myMass;
    Value myMask;
    ImageFct myMeasure;
    ImageFct myDistance2;
    Value myR2Max;
};

#endif
