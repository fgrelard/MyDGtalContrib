#ifndef EDGE_H
#define EDGE_H


#include <vector>
#include <algorithm>
#include "shapes/Curve.h"
#include "geometry/Distance.h"
#include "DGtal/topology/DomainMetricAdjacency.h"

template <typename Container>
class Edge : public Curve<Container> {
public:
    typedef typename Container::value_type Point;
public:
    Edge(const Container& aCurve) : Curve<Container>(aCurve) {}
    Edge(const Edge& other) : Curve<Container>(other) {}
public:
    std::vector< Edge<Container>* > neighboringEdges(const std::vector< Edge<Container>* >& edges,
                                                     const Container& branchingPoints);

};

template <typename Container>
std::vector< Edge<Container>* > Edge<Container>::neighboringEdges(const std::vector< Edge<Container>* >& edges,
                                                                  const Container& branchingPoints) {
    typedef typename Container::value_type Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::MetricAdjacency<Space, 3> MetricAdjacency;

    std::vector< Edge<Container>* > neighbors;
    Point branchPoint;
    for (const Point& b : branchingPoints) {
        if (std::find(this->pointSet().begin(), this->pointSet().end(), b) != this->pointSet().end() ) {
            branchPoint = b;
        }
    }

    std::vector<Point> nb;
    std::back_insert_iterator<std::vector<Point> > inserter(nb);
    MetricAdjacency::writeNeighbors(inserter, branchPoint);

    for (Edge<Container>* edge : edges) {
        Container setEdge = edge->pointSet();
        if (Distance::sameContainer(setEdge, this->pointSet())) continue;
        for (const Point& n : nb) {
            if (find(setEdge.begin(), setEdge.end(), n) != setEdge.end())
                neighbors.push_back(edge);
        }
    }

    return neighbors;
}

#endif
