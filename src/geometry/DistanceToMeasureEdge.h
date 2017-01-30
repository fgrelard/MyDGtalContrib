#ifndef DISTANCE_TO_MEASURE_EDGE_H
#define DISTANCE_TO_MEASURE_EDGE_H

#include "geometry/DistanceToMeasure.h"
#include "geometry/SetProcessor.h"
#include "Statistics.h"

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
                bool found = false;
                for (int i = 1; i < max; i++) {
                        m = DGtal::NumberTraits<Value>::ZERO;
                        d2 = DGtal::NumberTraits<Value>::ZERO;
                        Ball<Point> ball(p, i);
                        ImageFct intersection = ball.intersection(this->myMeasure);
                        DigitalSet empty(intersection.domain());
                        RealPoint g = Statistics<DigitalSet>(empty).centerOfMass(intersection);
                        RealVector dir = (g - p).getNormalized();
                        DigitalSet dirPlus = ball.pointsInHalfBall(dir);
                        DigitalSet dirMinus = ball.pointsInHalfBall(-dir);
                        for (const Point& pplus : dirPlus) {
                                if (!this->domain().isInside(pplus)) continue;
                                Value mpt = this->myMeasure( pplus );
                                double dist = l2Metric(p, pplus);
                                m += mpt;
                                d2 += mpt * dist * dist;
                        }
                        for (const Point& pminus : dirMinus) {
                                if (!this->domain().isInside(pminus)) continue;
                                Value mpt = this->myMeasure( pminus );
                                double dist = l2Metric(p, pminus);
                                m -= mpt;
                                d2 -= mpt * dist * dist;
                        }
                        if (m >= this->myMass) {
                                found = true;
                                break;
                        }

                }
                if (found)
                        return d2/m;
                return this->myR2Max;
        }

};

#endif
