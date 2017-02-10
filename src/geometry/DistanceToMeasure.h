#ifndef DISTANCE_TO_MEASURE_H
#define DISTANCE_TO_MEASURE_H

#include "DGtal/base/Trace.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "geometry/DistanceToPointFunctor.h"


template <typename ImageFct>
class DistanceToMeasure {
public:
        typedef typename ImageFct::Value   Value;
        typedef typename ImageFct::Point   Point;
        typedef typename ImageFct::Domain  Domain;
        typedef typename Domain::Space     Space;
        typedef typename Space::RealVector RealVector;
        typedef DGtal::MetricAdjacency<Space,1> Adjacency;

public:

        DistanceToMeasure( Value m0, const ImageFct& measure, Value rmax = 10.0, bool initialize = true )
                : myMass( m0 ), myMeasure( measure ), myDistance2( myMeasure.domain() ),
                  myR2Max( rmax*rmax ) {
                if (initialize)
                        init( );
        }


        DistanceToMeasure(const DistanceToMeasure& other) : myMass(other.myMass), myMeasure(other.myMeasure), myDistance2(other.myDistance2), myR2Max(other.myR2Max) {
        }

        virtual void init( ) {
                double       nb = myDistance2.domain().size();
                unsigned int i  = 0;

                for ( typename Domain::ConstIterator it = myDistance2.domain().begin(),
                              itE = myDistance2.domain().end(); it != itE; ++it, ++i )
                {
                        DGtal::trace.progressBar( i, nb );
                        myDistance2.setValue( *it, computeDistance2( *it ) );
                }
        }

        inline const ImageFct& measure() const {
                return myMeasure;
        }

        inline Value operator()( const Point& p ) const {
                return distance( p );
        }

        Value distance( const Point& p ) const {
                return sqrt( distance2( p ) );
        }


        Value distance2( const Point& p ) const {
                return myDistance2( p );
        }

        Value safeDistance2( const Point& p ) const {
                if ( myDistance2.domain().isInside( p ) )
                        return myDistance2( p );
                else return myDistance2( box( p ) );
        }

        Point box( const Point& p ) const {
                Point q = p.sup( myDistance2.domain().lowerBound() );
                return q.inf( myDistance2.domain().upperBound() );
        }

        virtual RealVector projection( const Point& p ) const {
                std::vector<Point> neighborsP;
                std::back_insert_iterator<std::vector<Point> > outIterator(neighborsP);
                Adjacency::writeNeighbors(outIterator, p);

                typedef typename std::vector<Point>::iterator Iterator;
                Value distance_center = distance2( p );
                RealVector vectorToReturn;
                for (Iterator it = neighborsP.begin(), ite = neighborsP.end();
                     it != ite; ++it) {
                        Point n = *it;
                        Value distance = (myDistance2.domain().isInside(*it)) ? distance2( n ) : distance_center;
                        for (int d = 0; d < Point::dimension; d++) {
                                if (p[d] != n[d]) {
                                        Point otherPoint = n;
                                        otherPoint[d] = p[d] + (p[d] - n[d]);
                                        Value otherDistance = (myDistance2.domain().isInside(otherPoint)) ? distance2( otherPoint ) : distance_center;
                                        if (otherPoint[d] > n[d]) {
                                                Value tmpDistance = otherDistance;
                                                otherDistance  = distance;
                                                distance = tmpDistance;
                                        }
                                        vectorToReturn[d] = ( std::fabs( distance - distance_center) >= std::fabs( distance_center - otherDistance) ) ? -(distance - distance_center) / 2.0 : -(distance_center - otherDistance) / 2.0;
                                }
                        }
                }
                return vectorToReturn;
        }

        virtual Value computeDistance2( const Point& p ) {
                typedef DGtal::ExactPredicateLpSeparableMetric<Space,2> Distance;
                typedef DistanceToPointFunctor<Distance>         DistanceToPoint;

                typedef DGtal::DistanceBreadthFirstVisitor< Adjacency, DistanceToPoint, std::set<Point> >
                        DistanceVisitor;
                typedef typename DistanceVisitor::Node MyNode;

                Value             m  = DGtal::NumberTraits<Value>::ZERO;
                Value             d2 = DGtal::NumberTraits<Value>::ZERO;
                Adjacency             graph;
                Distance l2;
                DistanceToPoint   d2pfct( l2, p );
                DistanceVisitor   visitor( graph, d2pfct, p );

                Value last = d2pfct( p );
                MyNode node;
                while ( ! visitor.finished() )
                {
                        node = visitor.current();
                        if ( ( node.second != last ) // all the vertices of the same layer have been processed.
                             && ( m >= myMass ) ) break;
                        if ( node.second > myR2Max ) { d2 = m * myR2Max; break; }
                        if ( myMeasure.domain().isInside( node.first ) )
                        {
                                Value mpt  = myMeasure( node.first );
                                d2        += mpt * node.second * node.second;
                                m         += mpt;
                                last       = node.second;
                                visitor.expand();
                        }
                        else
                                visitor.ignore();
                }
                return d2 / m;
        }

        Domain domain() const { return myMeasure.domain(); }

protected:
        Value myMass;
        const ImageFct& myMeasure;
        ImageFct myDistance2;
        Value myR2Max;
};

#endif
