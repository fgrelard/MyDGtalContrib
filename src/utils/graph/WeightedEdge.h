#ifndef WEIGHTED_EDGE_H
#define WEIGHTED_EDGE_H

#include "graph/Edge.h"

template <typename Container>
class WeightedEdge : public Edge<Container> {
public:
	typedef typename Container::value_type Point;
public:
    WeightedEdge(const Container& points, int label) : Edge<Container>(points), myLabel(label) {}
    WeightedEdge(const WeightedEdge& e) : Edge<Container>(e), myLabel(e.myLabel) {}
	int getLabel() const { return myLabel; }
	void setLabel(int label) { myLabel = label; }
    std::vector< WeightedEdge<Container>* > neighboringEdges(const std::vector< WeightedEdge<Container>* >& edges,
														  const std::vector<Point>& branchingPoints);
private:
	int myLabel;
};


template <typename Container>
std::vector< WeightedEdge<Container>* > WeightedEdge<Container>::neighboringEdges(const std::vector< WeightedEdge<Container>* >& edges,
																			const std::vector<Point>& branchingPoints) {
    typedef typename Container::value_type Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::MetricAdjacency<Space, 3> MetricAdjacency;

    std::vector< WeightedEdge<Container>* > neighbors;
    Point branchPoint;
    for (const Point& b : branchingPoints) {
        if (std::find(this->pointSet().begin(), this->pointSet().end(), b) != this->pointSet().end() ) {
            branchPoint = b;
        }
    }

    std::vector<Point> nb;
    std::back_insert_iterator<std::vector<Point> > inserter(nb);
    MetricAdjacency::writeNeighbors(inserter, branchPoint);

    for (WeightedEdge<Container>* edge : edges) {
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
