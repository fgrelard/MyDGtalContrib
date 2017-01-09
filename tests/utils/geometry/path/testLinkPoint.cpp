#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "geometry/path/LinkPointAlgorithm.h"
#include "geometry/path/BezierCasteljauLinkAlgorithm.h"
#include "geometry/path/BezierLinkAlgorithm.h"
#include "geometry/path/BresenhamAlgorithm.h"
#include "geometry/path/DijkstraAlgorithm.h"
#include "geometry/path/AStarAlgorithm.h"
#include "shapes/Curve.h"
#include "shapes/Border.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

using namespace DGtal;
using namespace std;

void testBezier(Viewer3D<>& viewer) {
    Z3i::Point first(2, 16, 31);
    Z3i::Point destination(3, -1, 22);
    Z3i::Point control2(0, 1, 25);
    Z3i::Point control1(1, 5, 10);

    BezierLinkAlgorithm<Z3i::Point> bezierAlgo(first, destination, control1, control2);

    Curve<vector<Z3i::Point> > points = bezierAlgo.linkPoints();
    viewer << CustomColors3D(Color::Yellow, Color::Yellow) << first << destination << control1 << control2;
    for (const Z3i::Point& p : points.pointSet()) {
        viewer << CustomColors3D(Color::Red, Color::Red) << p;
    }

    viewer << Viewer3D<>::updateDisplay;
}

void testBezierDeCasteljau(Viewer3D<>& viewer) {
    Z3i::Point first(2, 16, 31);
    Z3i::Point destination(3, -1, 23);
    Z3i::Point control1(0, 1, 25);
    Z3i::Point control2(1, 5, 25);

    BezierCasteljauLinkAlgorithm<Z3i::Point> bezierAlgo(first, destination, control1, control2);

    Curve<vector<Z3i::Point> > points = bezierAlgo.linkPoints();
    viewer << CustomColors3D(Color::Yellow, Color::Yellow) << first << destination << control1 << control2;
    for (const Z3i::Point& p : points.pointSet()) {
        viewer << CustomColors3D(Color::Red, Color::Red) << p;
    }

    viewer << Viewer3D<>::updateDisplay;
}

void testLinking() {
    using namespace Z3i;
    Point first(233, 276, 172);
    Point second(226, 276, 185);

    BresenhamAlgorithm<Z3i::Point> bresenhamAlgo(first, second);

    Curve<vector<Z3i::Point> > points = bresenhamAlgo.linkPoints();
    for (const Z3i::Point& p : points.pointSet()) {
        trace.info() << p << endl;
    }
}

void testDijkstra(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef Z3i::Object26_6 Graph;
        typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
        typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;

        Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        Border<Z3i::DigitalSet> border(setVolume);

        Z3i::DigitalSet setSurface = border.pointSet();
        Z3i::Point source(1, -9, 2);
        Z3i::Point destination(-9, -1, 8);
        DijkstraAlgorithm<Z3i::Point, Z3i::DigitalSet> dijAlgo(source, destination, setSurface);
        Curve<vector<Z3i::Point> > path = dijAlgo.linkPoints();


        const Color CURVE3D_COLOR( 100, 100, 140, 128 );

        for (const Z3i::Point& p : path.pointSet()) {
                viewer << CustomColors3D(Color::Green, Color::Green) << p;
        }
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << source << destination;
        viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << setSurface;
        viewer << Viewer3D<>::updateDisplay;
}

void testAStar(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef Z3i::Object26_6 Graph;
        typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
        typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;

        Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        Border<Z3i::DigitalSet> border(setVolume);

        Z3i::DigitalSet setSurface = border.pointSet();

        Z3i::Point source(1, -9, 2);
        Z3i::Point destination(-9, -1, 8);
        AStarAlgorithm<Z3i::Point, Z3i::DigitalSet> aAlgo(source, destination, setSurface);
        Curve<vector<Z3i::Point> > path = aAlgo.linkPoints();

        const Color CURVE3D_COLOR( 100, 100, 140, 128 );
        for (const Z3i::Point& p : path.pointSet()) {
                viewer << CustomColors3D(Color::Green, Color::Green) << p;
        }
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << source << destination;


        viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << setSurface;
        viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    //testBezier(viewer);
    //testBezierDeCasteljau(viewer);
//    testLinking();
//    testDijkstra(viewer);
    testAStar(viewer);
    app.exec();

    return 0;
}
