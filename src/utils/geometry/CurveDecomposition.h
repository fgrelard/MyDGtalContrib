#ifndef CURVE_DECOMPOSITION_H
#define CURVE_DECOMPOSITION_H

#include <vector>
#include <queue>
#include "graph/WeightedEdge.h"
#include "geometry/Curve.h"
#include "geometry/PointUtil.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"

template <typename Container>
class CurveDecomposition {
public:
    typedef typename Container::value_type Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;
    typedef typename Space::RealVector RealVector;
    typedef DGtal::HyperRectDomain< Space > Domain;
    typedef typename DGtal::DigitalSetSelector< Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS >::Type DigitalSet;

    typedef Curve<std::vector<Point> > CurveOrdered;
    typedef Edge<DigitalSet> GraphEdge;
    typedef WeightedEdge<DigitalSet> WeightedGraphEdge;


public:
    CurveDecomposition() = delete;
    CurveDecomposition(const Container& aCurve,
                       const Container& aBranchingPoints) : myCurve(aCurve),
                                                             myBranchingPoints(aBranchingPoints) {}
public:

    CurveOrdered curveTraversalForGraphDecomposition(const Point& start);

    std::vector< GraphEdge > constructGraph(const CurveOrdered& curve);

    std::vector< WeightedGraphEdge* > hierarchicalDecomposition(const std::vector< GraphEdge >& edges,
                                                                const std::vector< Point >& endPoints);

private:
    Curve<Container> myCurve;
    Container myBranchingPoints;
};

template <typename Container>
Curve< std::vector< typename CurveDecomposition<Container>::Point> >
CurveDecomposition<Container>::curveTraversalForGraphDecomposition(const typename CurveDecomposition<Container>::Point& p) {

    typedef DGtal::DepthFirstVisitor<ObjectType, std::set<Point> > Visitor;
    typedef typename Visitor::Node MyNode;
    typedef DGtal::ExactPredicateLpSeparableMetric<Space,2> L2Metric;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
    L2Metric l2Metric;
    std::vector<Point> existingSkeletonOrdered;
    ObjectType graph(dt26_6, myCurve.pointSet());
    Visitor visitor( graph, p );
    MyNode node;

    std::pair<Point, double> previous;
    while ( !visitor.finished() )
    {
        node = visitor.current();
        if (node.second != 0 && ((int)node.second - previous.second) <= 0) {
            std::vector<Point> neighbors;
            std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
            graph.writeNeighbors(inserter, node.first);
            double minDistance = std::numeric_limits<double>::max();
            Point cand;
            for (const Point& n : neighbors) {
                if (find(existingSkeletonOrdered.begin(), existingSkeletonOrdered.end(), n) != existingSkeletonOrdered.end()) {
                    double currentDistance = l2Metric(n, node.first);
                    if (currentDistance < minDistance) {
                        minDistance = currentDistance;
                        cand = n;
                    }
                }
            }
            existingSkeletonOrdered.push_back(cand);

        }
        previous = node;
        existingSkeletonOrdered.push_back(node.first);
        visitor.expand();
    }
    return existingSkeletonOrdered;
}

template <typename Container>
std::vector< typename CurveDecomposition<Container>::GraphEdge >
CurveDecomposition<Container>::
constructGraph(const Curve<std::vector<typename CurveDecomposition<Container>::Point> >& orderedCurve) {

    typedef DGtal::ExactPredicateLpSeparableMetric<Space,2> L2Metric;
    L2Metric l2Metric;

    std::vector< GraphEdge > graph;
    int index = 0;
    Point previous;

    Domain domain = PointUtil::computeBoundingBox<Domain>(orderedCurve.pointSet());
    DigitalSet toAdd(domain);
    std::vector<Point> orderedCurveSet = orderedCurve.pointSet();
    for (int i = 0, end = orderedCurveSet.size(); i < end; i++) {
        Point current = orderedCurveSet[i];
        if (l2Metric(previous,current) <= sqrt(3) || previous == Point())
            toAdd.insert(current);
        if (find(myBranchingPoints.begin(), myBranchingPoints.end(), current) != myBranchingPoints.end() ||
            l2Metric(previous,current) > sqrt(3)) {
            graph.push_back(GraphEdge(toAdd));
            index++;
            toAdd.clear();
            toAdd.insert(current);
        }
        previous = current;
    }
    return graph;
}


template <typename Container>
std::vector<typename CurveDecomposition<Container>::WeightedGraphEdge* >
CurveDecomposition<Container>::
hierarchicalDecomposition(const std::vector<typename CurveDecomposition<Container>::GraphEdge >& edges,
                          const std::vector<typename CurveDecomposition<Container>::Point>& endPoints) {

    std::queue<WeightedGraphEdge*> edgeQueue;
    std::vector<WeightedGraphEdge*> hierarchyGraph;
    for (const DigitalSet& edge : edges.pointSet()) {
        WeightedGraphEdge* levelEdge = new WeightedGraphEdge(edge, std::numeric_limits<int>::max());

        for (const Point& e : endPoints) {
            if (edge.find(e) != edge.end()) {
                levelEdge->setLabel(1);
                edgeQueue.push(levelEdge);
                break;
            }
        }
        hierarchyGraph.push_back(levelEdge);
    }

    while (!edgeQueue.empty()) {
        WeightedGraphEdge* edgeCurrent  = edgeQueue.front();
        DigitalSet edgeCurrentSet = edgeCurrent->pointSet();
        edgeQueue.pop();
        std::vector<WeightedGraphEdge*> neighborEdges = neighboringEdges(hierarchyGraph, edgeCurrentSet, myBranchingPoints);
        for (WeightedGraphEdge* neighCurrent : neighborEdges) {
            int label = edgeCurrent->getLabel()+1;
            if (neighCurrent->getLabel() > label) {
                neighCurrent->setLabel(label);
                edgeQueue.push(neighCurrent);
            }
        }
    }

    return hierarchyGraph;

}








#endif
