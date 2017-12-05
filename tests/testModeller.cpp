#include "Modeller.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/viewers/Viewer3D.h"

using namespace DGtal;
using namespace Z3i;

void testDrawBallVisitor(Viewer3D<>& viewer) {
    typedef DGtal::MetricAdjacency<Space, 1> Adjacency;

    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 1> Distance;
    typedef DistanceToPointFunctor<Distance> DistanceToPoint;

    typedef DGtal::DistanceBreadthFirstVisitor<Adjacency, DistanceToPoint, std::set<Point> >
        DistanceVisitor;
    typedef typename DistanceVisitor::Node MyNode;
    Adjacency graph;
    Distance l2;
    Point p;
    DistanceToPoint d2pfct(l2, p);
    DistanceVisitor visitor(graph, d2pfct, p);

    MyNode node;
    std::vector<MyNode> vec;
    visitor.current();
    visitor.getCurrentLayer(vec);
    visitor.expandLayer();
    visitor.getCurrentLayer(vec);
    visitor.expandLayer();
    visitor.getCurrentLayer(vec);
    visitor.expandLayer();
    visitor.getCurrentLayer(vec);
    visitor.expandLayer();
    visitor.getCurrentLayer(vec);
    node = visitor.current();
    DGtal::trace.info() << node.second << " " << vec.size() << std::endl;
    for (const auto& p : vec) {
        viewer << CustomColors3D(Color::Red, Color::Red) << p.first;
    }
}

DigitalSet translate(const DigitalSet& initial, const Point& translationVector) {
        DigitalSet out(Domain(initial.domain().lowerBound() + translationVector,
                              initial.domain().upperBound() + translationVector));
        for (const Point& p : initial) {
                out.insert(p+translationVector);
        }
        return out;
}

void testDrawCircle(Viewer3D<>& viewer) {
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet circle = modeller.drawCircle(5, Point(0, 0, 0));
        viewer << CustomColors3D(Color::Red, Color::Red) << circle;
}

void testDrawDisk(Viewer3D<>& viewer){
         Modeller<Z3i::DigitalSet> modeller;
         Z3i::DigitalSet circle = modeller.drawDisk(5, Point(11, 0, 0));
         viewer << CustomColors3D(Color::Red, Color::Red) << circle;
}


void testDrawCone(Viewer3D<>& viewer){
         Modeller<Z3i::DigitalSet> modeller;
         Z3i::DigitalSet cone = modeller.drawCone(Point(18, 0,0), RealVector(1,0,0), 5, 5);
         viewer << CustomColors3D(Color::Red, Color::Red) << cone;
}

void testDrawCylinder(Viewer3D<>& viewer){
         Modeller<Z3i::DigitalSet> modeller;
         Z3i::DigitalSet cylinder = modeller.drawCylinder(20, 5);
         Z3i::DigitalSet c = translate(cylinder, Point(30,0,0));
         viewer << CustomColors3D(Color::Red, Color::Red) << c;
}

void testDrawDeformedCylinder(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet cylinder = modeller.drawDeformedCylinder(20, 5);
        Z3i::DigitalSet c = translate(cylinder, Point(45,0,0));
        viewer << CustomColors3D(Color::Red, Color::Red) << c;

}

void testCreateHelixCurve(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet curve = modeller.createHelixCurve(20, 5, 3);
        Z3i::DigitalSet c = translate(curve, Point(60,0,0));
        viewer << CustomColors3D(Color::Red, Color::Red) << c;
}

void testCreateStraightLine(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet curve = modeller.createStraightLine(20);
        Z3i::DigitalSet c = translate(curve, Point(75,0,0));
        viewer << CustomColors3D(Color::Red, Color::Red) << c;
}
void testCreateSyntheticAirwayTree(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller(1);
         Z3i::DigitalSet curve(Z3i::Domain(Point(-100,-100,-100), Point(300,300,300)));
         modeller.createSyntheticAirwayTree(curve, 3, 20, 0, 0, Point(95,0,0));
//         Z3i::DigitalSet c = translate(curve, Point(95,0,0));
         viewer << CustomColors3D(Color::Red, Color::Red) << curve;
}

void testCreateLogarithmicCurve(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet curve = modeller.createLogarithmicCurve(20);
        Z3i::DigitalSet c = translate(curve, Point(105,0,0));
        viewer << CustomColors3D(Color::Red, Color::Red) << c;
}


void testCreateVolumeFromCurve(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet curve = modeller.createLogarithmicCurve(20);
        Z3i::DigitalSet c = translate(curve, Point(120,0,0));
        Z3i::DigitalSet vol = modeller.createVolumeFromCurve(c, 5);
        viewer << CustomColors3D(Color::Red, Color::Red) << vol;
}

void testCreateRotatedVolumeFromCurve(Viewer3D<>& viewer){
         Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet curve = modeller.createLogarithmicCurve(20);
        Z3i::DigitalSet c = translate(curve, Point(140,0,0));
        Z3i::DigitalSet vol = modeller.createRotatedVolumeFromCurve(c, 5, M_PI/4);
        viewer << CustomColors3D(Color::Red, Color::Red) << vol;
}

void testAddNoise(Viewer3D<>& viewer){
        Modeller<Z3i::DigitalSet> modeller;
        Z3i::DigitalSet cylinder = modeller.drawCylinder(20, 5);
        Z3i::DigitalSet curve = modeller.addNoise(cylinder, 0.3);
        Z3i::DigitalSet c = translate(curve, Point(160,0,0));
        viewer << c;
}

void testDrawSurface(Viewer3D<>& viewer) {
    Modeller<Z3i::DigitalSet> modeller;
    Z3i::DigitalSet surface = modeller.drawSurface(100);
    viewer << surface;
}

void testCreateDeformedSyntheticAirwayTree(Viewer3D<>& viewer) {
    Modeller<Z3i::DigitalSet> modeller;
    Z3i::DigitalSet airway(Z3i::Domain(Z3i::Point(0,0,0), Z3i::Point(100,100,100)));
    Z3i::Point p(0,0,0);
    modeller.createDeformedSyntheticAirwayTree(airway, 4, 50, 0, 0, p);
    viewer << airway;
}

int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testCreateDeformedSyntheticAirwayTree(viewer);
//        testDrawSurface(viewer);
//        testDrawBallVisitor(viewer);
        // testDrawCircle(viewer);
        // testDrawDisk(viewer);
        // testDrawCone(viewer);
        // testDrawCylinder(viewer);
        // testDrawDeformedCylinder(viewer);
        // testCreateHelixCurve(viewer);
        // testCreateStraightLine(viewer);
        // testCreateSyntheticAirwayTree(viewer);
        // testCreateLogarithmicCurve(viewer);
        // testCreateVolumeFromCurve(viewer);
        // testCreateRotatedVolumeFromCurve(viewer);
        // testAddNoise(viewer);
        viewer << Viewer3D<>::updateDisplay;
        app.exec();
}
