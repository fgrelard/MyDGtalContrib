#ifndef RECENTER_SKELETON_H
#define RECENTER_SKELETON_H

#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "vcm/CuttingPlaneEstimator.h"
#include "geometry/CurveDecomposition.h"

template <typename Container>
class RecenterSkeleton {
public:
        typedef typename Container::Point Point;
        typedef typename Container::Space Space;
        typedef typename Space::RealVector RealVector;
        typedef typename DGtal::functors::BallConstantPointFunction<Point, double> KernelFunction;
        typedef typename DGtal::HyperRectDomain< Space > Domain;
        typedef OrthogonalPlaneEstimator<Container, KernelFunction> PlaneEstimator;
        typedef typename PlaneEstimator::Plane Plane;
        typedef typename PlaneEstimator::L2Metric L2Metric;
        typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;
        typedef Edge<Container> GraphEdge;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;
        typedef CuttingPlaneEstimator<PlaneEstimator> CuttingPlane;

public:
        RecenterSkeleton() = delete;
        RecenterSkeleton(const Container& skeleton,
                         const Container& volume,
                         const DTL2& aDT);
        ~RecenterSkeleton();
        RecenterSkeleton(const RecenterSkeleton& other);

public:
        Container recenter();
public:
        std::vector<Container> adjacentBranchesToPoint(const Point& p,
                                                       const std::vector<GraphEdge>& graph);

        std::vector<Container> restrictBranches(const std::vector<Container>& branches, const Point& p);
        Container referenceBranch(const std::vector<Container>& branches, const std::vector<Container>& subSetBranches);
        Container planeToBranch(const Plane& plane, const std::vector<Container>& branches);
        std::vector<Plane> alignPlanes();
        Container subVolume();

private:
        Container* mySkeleton;
        Container* myVolume;
        DTL2* myDT;

        //Internals
private:
        PlaneEstimator* myPlaneEstimator;

};

template <typename Container>
RecenterSkeleton<Container>::
RecenterSkeleton(const Container& skeleton,
                 const Container& volume,
                 const DTL2& aDT) {
        mySkeleton = new Container ( skeleton );
        myVolume = new Container ( volume );
        myDT = new DTL2 ( aDT );
        double r = 5;
        double R = 30;
        KernelFunction chi(1.0, r);
        int connexity = 6;
        myPlaneEstimator = new PlaneEstimator(volume, chi, R, r, connexity );
}

template <typename Container>
RecenterSkeleton<Container>::
~RecenterSkeleton() {
        if (mySkeleton) {
                delete mySkeleton;
                mySkeleton = 0;
        }
        if (myVolume) {
                delete myVolume;
                myVolume = 0;
        }
        if (myDT) {
                delete myDT;
                myDT = 0;
        }
        if (myPlaneEstimator) {
                delete myPlaneEstimator;
                myPlaneEstimator = 0;
        }
}

template <typename Container>
RecenterSkeleton<Container>::
RecenterSkeleton(const RecenterSkeleton& other) {
        mySkeleton = new Container ( *other.mySkeleton );
        myVolume = new Container ( *other.myVolume );
        myPlaneEstimator = new PlaneEstimator( *other.myPlaneEstimator );
        myDT = new DTL2( *other.myDT );
}

template <typename Container>
Container
RecenterSkeleton<Container>::
recenter() {
        Container branching = CurveProcessor<Container>(*mySkeleton).branchingPoints();
        std::vector<GraphEdge> graph = CurveDecomposition<Container>(*mySkeleton, branching).branchDecomposition();
        for (const Point& b : branching) {
                std::vector<Container> adjacentEdges = adjacentBranchesToPoint(b, graph);
                double radius = max_element(adjacentEdges.begin(), adjacentEdges.end(), [&](const Container& e1,
                                                                                            const Container& e2) {
                                                    return e1.size() < e2.size();
                                            })->size()*0.5;

                std::vector<Container> restricted = restrictBranches(adjacentEdges, b);
                CuttingPlane cuttingPE(restricted, *myPlaneEstimator, *myVolume, *myDT);
                std::vector<Plane> cuttingPlanes = cuttingPE.cuttingPlanes();
                if (cuttingPlanes.size() < 2) continue;
        }

}


template <typename Container>
std::vector<Container>
RecenterSkeleton<Container>::
adjacentBranchesToPoint(const Point& p,
                     const std::vector<GraphEdge>& graph) {
        std::vector<Container> adjacent;
        Adj26 adj26;
        Adj6 adj6;
        DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
        ObjectType objSkel(dt26_6, *mySkeleton);
        std::vector<Point> neigh;
        std::back_insert_iterator<std::vector<Point>> inserter(neigh);
        objSkel.writeNeighbors(inserter, p);

        for (const Container& edge : graph) {
                for (const Point& n : neigh) {
                        if (edge.find(n) != edge.end()) {
                                adjacent.push_back(edge);
                                break;
                        }
                }
        }
        return adjacent;
}

template <typename Container>
Container
RecenterSkeleton<Container>::
planeToBranch(const Plane& plane,
              const std::vector<Container>& branches) {
        auto iterator = branches.find(plane.getCenter());
        if (iterator != branches.end())
                return *iterator;
        return Container(Domain(Point::zero, Point::zero));
}
template <typename Container>
Container
RecenterSkeleton<Container>::
referenceBranch(const std::vector<Container>& branches, const std::vector<Container>& subSetBranches) {
        Container restrictEdgeB(myVolume->domain());
        for (const Container& branch : branches) {
                unsigned int cpt = 0;
                for (const Container& branchSub : subSetBranches) {
                        if (!SetProcessor<Container>(branch).sameContainer(branchSub)) {
                                cpt++;
                        }
                }
                if (cpt == subSetBranches.size())
                        restrictEdgeB = branch;
        }
        return restrictEdgeB;
}

template <typename Container>
std::vector<Container>
RecenterSkeleton<Container>::
restrictBranches(const std::vector<Container>& branches,
                 const Point& p) {
        std::vector<Container> restrictedBranches;
        for (const Container& branch : branches) {
                Container restricted = CurveProcessor<Container>(branch).subCurve(p, branch.size());
                restrictedBranches.push_back(restricted);
        }
        return restrictedBranches;
}


#endif
