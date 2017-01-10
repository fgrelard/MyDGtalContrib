#include "graph/Edge.h"
#include "graph/WeightedEdge.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"

using namespace DGtal;
using namespace std;


void testEdge() {
        trace.beginBlock("Test edge");
        Z3i::DigitalSet aSet(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                         Z3i::Point(100,100,100)));
        aSet.insert(Z3i::Point(0,0,0));
        aSet.insert(Z3i::Point(1,1,1));
        aSet.insert(Z3i::Point(1,1,2));
        Edge<Z3i::DigitalSet> curve(aSet);
        for (const Z3i::Point& p : curve)
                trace.info() << p << endl;
        trace.endBlock();
}

void testNeighborEdge() {
        trace.beginBlock("Test neighbor edge");
        Z3i::DigitalSet aSet(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                         Z3i::Point(100,100,100)));
        aSet.insert(Z3i::Point(0,0,0));
        aSet.insert(Z3i::Point(1,1,1));
        aSet.insert(Z3i::Point(1,1,2));
        Edge<Z3i::DigitalSet> curve(aSet);

        Z3i::DigitalSet aSet2(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                         Z3i::Point(100,100,100)));
        aSet2.insert(Z3i::Point(0,0,0));
        aSet2.insert(Z3i::Point(-1,-1,-1));
        aSet2.insert(Z3i::Point(-1,-1,-2));

        Z3i::DigitalSet branchingPoints(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                                     Z3i::Point(100,100,100)));

        branchingPoints.insert(Z3i::Point(0,0,0));
        Edge<Z3i::DigitalSet> curve2(aSet2);
        std::vector<Edge<Z3i::DigitalSet>* > edges;
        edges.push_back(new Edge<Z3i::DigitalSet>(curve));
        edges.push_back(new Edge<Z3i::DigitalSet>(curve2));

        auto v = curve.neighboringEdges(edges, branchingPoints);
        trace.info() << v.size() << endl;
        trace.endBlock();
}

void testNeighborGraphEdge() {
        trace.beginBlock("Test neighbor graph edge");
        Z3i::DigitalSet aSet(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                         Z3i::Point(100,100,100)));
        aSet.insert(Z3i::Point(0,0,0));
        aSet.insert(Z3i::Point(1,1,1));
        aSet.insert(Z3i::Point(1,1,2));
        WeightedEdge<Z3i::DigitalSet> curve(aSet, 1);

        Z3i::DigitalSet aSet2(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                         Z3i::Point(100,100,100)));
        aSet2.insert(Z3i::Point(0,0,0));
        aSet2.insert(Z3i::Point(-1,-1,-1));
        aSet2.insert(Z3i::Point(-1,-1,-2));

        Z3i::DigitalSet branchingPoints(Z3i::Domain(Z3i::Point(-100,-100,-100),
                                                     Z3i::Point(100,100,100)));

        branchingPoints.insert(Z3i::Point(0,0,0));
        WeightedEdge<Z3i::DigitalSet> curve2(aSet2, 2);
        std::vector<WeightedEdge<Z3i::DigitalSet>* > edges;
        edges.push_back(new WeightedEdge<Z3i::DigitalSet>(curve));
        edges.push_back(new WeightedEdge<Z3i::DigitalSet>(curve2));

        auto v = curve.neighboringEdges(edges, branchingPoints);
        trace.info() << v.size() << endl;
        trace.endBlock();
}


int main() {
        testEdge();
        testNeighborEdge();
        testNeighborGraphEdge();
        return 0;
}
