#ifndef CURVE_PROCESSOR_H
#define CURVE_PROCESSOR_H

#include <vector>
#include <queue>

#include "shapes/Curve.h"
#include "DGtal/base/Common.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"

template <typename Container>
class CurveProcessor {

	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

public:
	typedef typename Container::value_type Point;
	typedef DGtal::MetricAdjacency<typename Container::Space, 1> Adj6;
	typedef DGtal::MetricAdjacency<typename Container::Space, 3> Adj26;
	typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
	typedef DGtal::Object<DT26_6, Container> ObjectType;
	typedef typename Container::Space::RealVector RealVector;

public:
	CurveProcessor(const Container& container, bool isOrdered = false) : myCurve(container) {}

public:

	Container ensureConnexity();

	Container endPoints();

	Container branchingPoints();

	Curve< std::vector<Point> > convertToOrderedCurve();

	Curve< std::vector<Point> > convertToOrderedCurve(const Point& startingPoint);



private:
    Curve<Container> myCurve;
};




template <typename Container>
Container CurveProcessor<Container>::ensureConnexity() {

	Adj26 adj26;
	Adj6 adj6;
	DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );
	Container cleanSet(myCurve.pointSet().domain());
	ObjectType obj(dt26_6, myCurve.pointSet());
	Container & S = obj.pointSet();
	cleanSet = S;
	for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
		ObjectType obj(dt26_6, cleanSet);
		if (obj.isSimple(*it)) {
		    cleanSet.erase(*it);
		}
	}

    return cleanSet;
}

template <typename Container>
Container CurveProcessor<Container>::endPoints() {

	Adj26 adj26;
	Adj6 adj6;
	DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );

	Container set = myCurve.pointSet();
    ObjectType objectSet(dt26_6, set);
    Container endPoints(set.domain());
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
	    Point p = *it;
		std::vector<Point> neighbors;
		std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
		objectSet.writeNeighbors(inserter, p);
		if (neighbors.size() <= 1)
			endPoints.insert(p);

		//Is it in same quadrant: case connectivity != 26
		else {
		    RealVector previous;
			bool isEndPoint = true;
			std::vector<RealVector> vectors;
			for (const Point& n : neighbors) {
			    RealVector dir = (n - p).getNormalized();
				vectors.push_back(dir);
			}
			//Min angle (resp max dot product) determined by two points with one varying coordinate
			for (int i = 0; i < vectors.size(); i++) {
				for (int j = i+1; j < vectors.size(); j++) {
					if (vectors[i].dot(vectors[j]) <= (1/(1+sqrt(2))) )
						isEndPoint = false;
				}
			}
			if (isEndPoint)
				endPoints.insert(p);
		}
	}
	return endPoints;
}



template <typename Container>
Container CurveProcessor<Container>::branchingPoints() {

	Adj26 adj26;
	Adj6 adj6;
	DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );

	Container set = myCurve.pointSet();

    ObjectType obj(dt26_6, set);
	Container criticalPoints(set.domain());
	for (const Point& s : set) {
		std::vector<Point> neighbors;
		std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
		obj.writeNeighbors(inserter, s);
		if (neighbors.size() > 2) {
			criticalPoints.insert(s);
		}
	}
	return criticalPoints;
}

template <typename Container>
Curve< std::vector< typename CurveProcessor<Container>::Point > > CurveProcessor<Container>::convertToOrderedCurve() {
	std::vector<Point> orientedEdge;

	Container set = myCurve.pointSet();
	if (set.size() == 0) return orientedEdge;

	Container e = endPoints();
    Point start = *(e.begin());
	return convertToOrderedCurve(start);

}


template <typename Container>
Curve< std::vector< typename CurveProcessor<Container>::Point > > CurveProcessor<Container>::convertToOrderedCurve(const typename CurveProcessor<Container>::Point& startingPoint) {
	std::vector<Point> orientedEdge;
	Container edge = myCurve.pointSet();
	if (edge.size() == 0) return orientedEdge;

	Adj26 adj26;
	Adj6 adj6;
	DT26_6 dt26_6 (adj26, adj6, DGtal::JORDAN_DT );

	ObjectType objEdge(dt26_6, edge);
	Point start = startingPoint;
	orientedEdge.push_back(start);
	bool toAdd = true;
	while (toAdd) {
		std::vector<Point> neighbors;
		std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
		objEdge.writeNeighbors(inserter, start);
		unsigned int cpt = 0;
		for (const Point& n : neighbors) {
			if (std::find(orientedEdge.begin(), orientedEdge.end(), n) == orientedEdge.end()) {
				orientedEdge.push_back(n);
				start = n;
			}
			else
				cpt++;
		}
		if (cpt == neighbors.size())
			toAdd = false;
	}
	return orientedEdge;
}





#endif
