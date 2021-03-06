#ifndef DISTANCE_TO_MEASURE_RELEVANT_SCALE_H
#define DISTANCE_TO_MEASURE_RELEVANT_SCALE_H

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
class DistanceToMeasureRelevantScale {
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

    DistanceToMeasureRelevantScale(Value m0, const ImageFct &measure, Value rmax = 10.0, Value mask = 1.1, bool initialize =
                                   true, bool speedup = true)
            : myMass(m0), myMeasure(measure), myDistance2(myMeasure.domain()),
              myR2Max(rmax * rmax), myMask(mask), mySpeedUp(speedup) {
        if (initialize)
            init();
    }

    DistanceToMeasureRelevantScale(const ImageFct& distance2) : myDistance2(distance2), myMeasure(distance2.domain()) {}


    DistanceToMeasureRelevantScale(const DistanceToMeasureRelevantScale &other) : myMass(other.myMass), myMeasure(other.myMeasure),
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
            if (myMeasure(current) != myMask)
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

    virtual Value maxAlongProj(const Point& p) const {
        RealVector proj = -projection(p);
        if (proj.norm() == 0) return 0;
        proj = proj.getNormalized();
        Value previousValue = std::sqrt(myDistance2(p));
        Value nextValue = previousValue;
        if (proj.norm() == 0 || previousValue == 0) return 0;

        Point current = p, previous = p;
        do {
            previousValue = nextValue;
            previous = current;
            int i = 1;
            while (current == previous) {
                current = previous + proj * i;
                i++;
            }
            if (!myDistance2.domain().isInside(current)) break;
            nextValue = std::sqrt(myDistance2(current));


        } while (previousValue <= nextValue);
        return previousValue;
    }

    virtual Value computeDistance2(const Point &p) {
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> Distance;
        typedef DistanceToPointFunctor<Distance> DistanceToPoint;

        typedef DGtal::DistanceBreadthFirstVisitor<Adjacency, DistanceToPoint, std::set<Point> >
                DistanceVisitor;
        typedef typename DistanceVisitor::Node MyNode;

        Value aMass = myMass * myMass;
        Value m = DGtal::NumberTraits<Value>::ZERO;
        Adjacency graph;
        Distance l2;
        DistanceToPoint d2pfct(l2, p);
        DistanceVisitor visitor(graph, d2pfct, p);
        Value last = d2pfct(p);
        MyNode node;
        Value firstMass = myMeasure(p);
        DGtal::Statistic<Value> stat(true);
        //stat.addValue(firstMass);

        Value previousMean = DGtal::NumberTraits<Value>::ZERO;
        Value currentMean = std::numeric_limits<Value>::min();
        while (!visitor.finished()) {
            node = visitor.current();
            std::vector<MyNode> vec;
            visitor.getCurrentLayer(vec);
            for (const MyNode &n : vec) {
                if (!myMeasure.domain().isInside(n.first)) continue;
                double currentColor = myMeasure(n.first);
                stat.addValue(currentColor);
            }
            if (mySpeedUp && node.second > 4 && stat.mean() >= previousMean) {
                currentMean = stat.mean();
                break;
            }
            m = stat.variance();
            if (node.second * node.second >= 3) {
                if (m >= aMass) {
                    currentMean = stat.mean();
                    break;
                }
                else if (node.second * node.second > myR2Max) {
                    currentMean = previousMean;
                    break;
                }
            } else if (m >= aMass) {
                break;
            }

            //Next layer
            previousMean = stat.mean();

            //Next layer
            visitor.expandLayer();
        }

        if (m == DGtal::NumberTraits<Value>::ZERO)
            return DGtal::NumberTraits<Value>::ZERO;
        if (currentMean >= previousMean)
            return DGtal::NumberTraits<Value>::ZERO;


        return node.second * node.second;
    }

    Domain domain() const { return myDistance2.domain(); }

protected:
    Value myMass;
    Value myMask;
    const ImageFct &myMeasure;
    ImageFct myDistance2;
    Value myR2Max;
    bool mySpeedUp = false;
};

#endif
