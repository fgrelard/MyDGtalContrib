#ifndef PRUNING_WITH_ORTHOGONAL_PLANES_H
#define PRUNING_WITH_ORTHOGONAL_PLANES_H

#include "vcm/OrthogonalPlaneEstimator.h"
#include "geometry/CurveDecomposition.h"
#include "geometry/CurveProcessor.h"
#include "graph/WeightedEdge.h"
#include "DGtal/kernel/CSpace.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"


template <typename Container>
class PruningWithOrthogonalPlanes {

        BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));
public:
        typedef typename Container::Point Point;
        typedef typename Container::Space Space;
        typedef typename Space::RealVector RealVector;
        typedef typename DGtal::functors::BallConstantPointFunction<Point, double> KernelFunction;
        typedef typename DGtal::HyperRectDomain< Space > Domain;
        typedef OrthogonalPlaneEstimator<Container, KernelFunction> PlaneEstimator;
        typedef typename PlaneEstimator::Plane Plane;
        typedef WeightedEdge<Container> GraphEdge;
        typedef typename PlaneEstimator::L2Metric L2Metric;
        typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;

public:
        PruningWithOrthogonalPlanes() = delete;
        PruningWithOrthogonalPlanes(const Container& skeleton,
                                    const Container& volume,
                                    double threshold);
        ~PruningWithOrthogonalPlanes();
        PruningWithOrthogonalPlanes(const PruningWithOrthogonalPlanes& other);

public:
        Container prune();

        double significanceMeasure(const Point& p);

        Container pruneEdgeTopologyPreserving(const Container& prunedSkeleton,
                                              double significanceMeasure,
                                              const GraphEdge& graphEdge);

private:
        Container* mySkeleton;
        Container* myVolume;

        DTL2* myDT;

        PlaneEstimator* myPlaneEstimatorCurve;
        PlaneEstimator* myPlaneEstimatorVol;

        Container myBranchingPoints;
        double myThreshold;
};

template <typename Container>
PruningWithOrthogonalPlanes<Container>::
PruningWithOrthogonalPlanes(const Container& skeleton,
                            const Container& volume, double threshold) : myBranchingPoints(Container(Domain(Point::zero,
                                                                                                            Point::zero))){
        mySkeleton = new Container ( skeleton );
        myVolume = new Container ( volume );
        double r = 5;
        double R = 30;
        KernelFunction chi(1.0, r);
        myPlaneEstimatorVol = new PlaneEstimator ( volume, chi, R, r );

        myPlaneEstimatorCurve = new PlaneEstimator( skeleton, chi, R, r);
        L2Metric l2Metric;
        myDT = new DTL2(volume.domain(), volume, l2Metric);

        CurveProcessor<Container> curveProc(*mySkeleton);
        myBranchingPoints = curveProc.branchingPoints();

        myThreshold = threshold * M_PI / 180.0 ;
}

template <typename Container>
PruningWithOrthogonalPlanes<Container>::
~PruningWithOrthogonalPlanes() {
        if (mySkeleton != 0) {
                delete mySkeleton;
                mySkeleton = 0;
        }
        if (myVolume != 0) {
                delete myVolume;
                myVolume = 0;
        }
        if (myPlaneEstimatorCurve != 0) {
                delete myPlaneEstimatorCurve;
                myPlaneEstimatorCurve = 0;
        }
        if (myPlaneEstimatorVol != 0) {
                delete myPlaneEstimatorVol;
                myPlaneEstimatorVol = 0;
        }
        if (myDT != 0) {
                delete myDT;
                myDT = 0;
        }

}

template <typename Container>
PruningWithOrthogonalPlanes<Container>::
PruningWithOrthogonalPlanes(const PruningWithOrthogonalPlanes& other) {
        mySkeleton = new Container ( *other.mySkeleton );
        myVolume = new Container ( *other.myVolume );
        myPlaneEstimatorCurve = new PlaneEstimator( *other.myPlaneEstimatorCurve );
        myPlaneEstimatorVol = new PlaneEstimator ( *other.myPlaneEstimatorVol );
        myDT = new DTL2( *other.myDT );
        myBranchingPoints = other.myBranchingPoints;
        myThreshold = other.myThreshold;

}

template <typename Container>
Container
PruningWithOrthogonalPlanes<Container>::
prune() {
        DGtal::trace.beginBlock("Pruning skeleton");
        CurveDecomposition<Container> curveDecompo(*mySkeleton, myBranchingPoints);
        std::vector<GraphEdge*> hierarchicalGraph = curveDecompo.graphDecomposition();
        int previousNumber = -1;
        Container prunedSkeleton = *mySkeleton;
        while (hierarchicalGraph.size() != previousNumber) {
                previousNumber = hierarchicalGraph.size();
                for (GraphEdge* graphEdge : hierarchicalGraph) {
                        if (graphEdge->size() == 0) continue;
                        double radius = graphEdge->size() * 0.4;
                        KernelFunction chi(1.0, radius);
                        myPlaneEstimatorCurve = new PlaneEstimator (*graphEdge, chi, 30, radius);
                        double sumAngle = 0;
                        std::for_each(graphEdge->begin(),
                                      graphEdge->end(),
                                      [&](const Point& p) {
                                              sumAngle += significanceMeasure(p);
                                      });
                        sumAngle /= graphEdge->size();
                        if (sumAngle > myThreshold) {
                                prunedSkeleton = pruneEdgeTopologyPreserving(prunedSkeleton, sumAngle, *graphEdge);

                        }
                }
                CurveProcessor<Container> curveProc(prunedSkeleton);
                myBranchingPoints = curveProc.branchingPoints();
                curveDecompo = CurveDecomposition<Container>(*mySkeleton, myBranchingPoints);
                hierarchicalGraph = curveDecompo.graphDecomposition();

        }
        return prunedSkeleton;
        DGtal::trace.endBlock();
}

template <typename Container>
double
PruningWithOrthogonalPlanes<Container>::
significanceMeasure(const Point& p) {

        Plane planeCurve = myPlaneEstimatorCurve->planeAt(p);

        double radiusVol = (*myDT)(p) + 2;
        myPlaneEstimatorVol->setRadius(radiusVol);
        Plane planeVol = myPlaneEstimatorCurve->convergentPlaneAt(p, *myVolume, radiusVol*10);
        RealVector normalCurve = planeCurve.getPlaneEquation().normal();
        RealVector normalVol = planeVol.getPlaneEquation().normal();
        double angle = normalVol.cosineSimilarity(normalCurve);
        double otherAngle = normalVol.cosineSimilarity(-normalCurve);
        angle = (angle < otherAngle) ? angle : otherAngle;
        return angle;
}

template <typename Container>
Container
PruningWithOrthogonalPlanes<Container>::
pruneEdgeTopologyPreserving(const Container& prunedSkeleton,
                            double significanceMeasure,
                            const GraphEdge& graphEdge) {
        Container difference(prunedSkeleton.domain());
        for (const Point& p : prunedSkeleton) {
                if (graphEdge.find(p) == graphEdge.end() ||
                    myBranchingPoints.find(p) != myBranchingPoints.end())
                        difference.insert(p);
        }
        Adj26 adj26;
        Adj6 adj6;
        DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
        ObjectType obj(dt26_6, difference);
        std::vector<ObjectType> objects;
        std::back_insert_iterator< std::vector<ObjectType> > inserter(objects);
        unsigned int nbCC = obj.writeComponents(inserter);
        if (nbCC == 1)
                return difference;
        return prunedSkeleton;

}


#endif
