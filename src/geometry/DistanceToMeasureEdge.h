#ifndef DISTANCE_TO_MEASURE_EDGE_H
#define DISTANCE_TO_MEASURE_EDGE_H

#include "geometry/DistanceToMeasure.h"
#include "geometry/SetProcessor.h"
#include "DGtal/math/Statistic.h"
#include "ShapeDescriptor.h"

template <typename ImageFct>
class DistanceToMeasureEdge : public DistanceToMeasure<ImageFct> {
public:
        typedef DistanceToMeasure<ImageFct> Base;
        typedef typename Base::Value Value;
        typedef typename Base::Point Point;
        typedef typename Base::Domain Domain;
        typedef typename Base::Space Space;
        typedef typename Space::RealPoint RealPoint;
        typedef typename Base::RealVector RealVector;
        typedef typename DGtal::DigitalSetSelector< Domain, DGtal::BIG_DS+DGtal::HIGH_BEL_DS >::Type DigitalSet;
        typedef DGtal::MetricAdjacency<Space, 2> Adjacency;

public:
        DistanceToMeasureEdge(Value m0, const ImageFct& measure,
                              Value rmax = 10.0) : Base(m0, measure, rmax, false) {

                init();
        }

public:
        void init( ) {
                double       nb = this->domain().size();
                unsigned int i  = 0;

                for ( typename Domain::ConstIterator it = this->domain().begin(),
                              itE = this->domain().end(); it != itE; ++it, ++i )
                {
                        DGtal::trace.progressBar( i, nb );
                        this->myDistance2.setValue( *it, computeDistance2( *it ) );
                }
        }

        Value computeDistance2( const Point& p) {
                typedef DGtal::ExactPredicateLpSeparableMetric<Space,2> Distance;

                Value             m  = DGtal::NumberTraits<Value>::ZERO;
                Value             d2 = DGtal::NumberTraits<Value>::ZERO;

                Distance l2Metric;
                double max = sqrt(this->myR2Max);
                int i = max;
                for (; i > 1; i--) {
                        m = DGtal::NumberTraits<Value>::ZERO;
                        d2 = DGtal::NumberTraits<Value>::ZERO;
                        Ball<Point> ball(p, i);
                        ImageFct intersection = ball.intersection(this->myMeasure);
                        DigitalSet empty(intersection.domain());
                        RealPoint g = ShapeDescriptor<DigitalSet>(empty).centerOfMass(intersection);
                        RealVector dir = (g - p).getNormalized();
                        DigitalSet dirPlus = ball.pointsInHalfBall(dir);
                        DigitalSet dirMinus = ball.pointsInHalfBall(-dir);
                        DigitalPlane<Space> plane(p, dir, 8);
                        DigitalSet pointsBall = ball.pointSet();
                        DigitalSet pointsInPlane = empty;//plane.intersectionWithSetOneCC(pointsBall);
                        std::vector<double> values;
                        for (const Point& pplus : dirPlus) {
                                if (!this->domain().isInside(pplus) ||
                                    pointsInPlane.find(pplus) != pointsInPlane.end()) continue;
                                Value mpt = this->myMeasure( pplus );
                                double dist = l2Metric(p, pplus);
                                m += mpt;
                                d2 += mpt * dist * dist;
                                values.push_back(m);
                        }
                        DGtal::Statistic<double> stats;
                        stats.addValues(values.begin(), values.end());
                        double std = sqrt(stats.variance());
                        m *= std;
                        for (const Point& pminus : dirMinus) {
                                if (!this->domain().isInside(pminus) ||
                                    pointsInPlane.find(pminus) != pointsInPlane.end()) continue;
                                Value mpt = this->myMeasure( pminus );
                                double dist = l2Metric(p, pminus);
                                m -= mpt;
                                d2 -= mpt * dist * dist;
                        }

                        if (m < this->myMass) {
                                break;
                        }

                }
                // m = std::fabs(m);
                // d2 = std::fabs(d2);
                if (m < 1)
                        m = 1;
                if (i != max)
                        return d2/m;
                return this->myR2Max;
        }

};

#endif
