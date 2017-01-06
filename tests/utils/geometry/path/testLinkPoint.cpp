#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "geometry/path/LinkPointAlgorithm.h"
#include "geometry/path/BezierCasteljauLinkAlgorithm.h"
#include "geometry/path/BezierLinkAlgorithm.h"
#include "geometry/path/BresenhamAlgorithm.h"
#include "geometry/Curve.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/io/viewers/Viewer3D.h"

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





int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    testBezier(viewer);
    //testBezierDeCasteljau(viewer);
    testLinking();
//	testDSSLinking();
    app.exec();

    return 0;
}
