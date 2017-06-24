#include <shapes/Polygon.h>
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

using namespace DGtal;
using namespace std;


void testPolygonInside(Viewer3D<>& viewer) {
//    Z2i::Point p1(0,0);
//    Z2i::Point p2(1,10);
//    Z2i::Point p3(0,-7);
//    Z2i::Point p4(1,-3);
    Z2i::Point p1(276, 189);
    Z2i::Point p2(290, 189);

    Polygon<Z2i::Point> polygon{p1, p2};
    Z2i::Domain d(Z2i::Point(275, 187), Z2i::Point(291, 191));
    for (const Z2i::Point& p : d) {
        Z3i::Point p3D(p[0], p[1], 0);
        if (polygon.isInside(p))
            viewer << CustomColors3D(Color::Red, Color::Red) << p3D;
    }
    viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testPolygonInside(viewer);
        app.exec();
        return 0;
}
