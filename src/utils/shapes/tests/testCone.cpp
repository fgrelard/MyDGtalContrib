#include "../Cone.h"

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include <QtGui/qapplication.h>

using namespace std;
using namespace DGtal;

void testCone(int argc, char** argv) {
            QApplication app(argc, argv);
            Viewer3D<> viewer;
            viewer.show();
            Z3i::Point center(0,0,0);
            Z3i::RealVector direction(-1,-1,-1);
            Cone<Z3i::Point, Z3i::RealVector> cone(center, direction.getNormalized(), 10, 100);
            bool contains1 = cone.contains(Z3i::Point(-2,-2,-2));
            bool contains2 = cone.contains(Z3i::Point(0,0,0));
            bool contains3 = cone.contains(Z3i::Point(0,0,1));
            trace.info() << contains1 << contains2 << contains3 << std::endl;
            std::vector<Z3i::Point> points = cone.pointSet();
            for (const auto& p : points) {
                    viewer << CustomColors3D(Color::Red, Color::Red) << p;
            }
            viewer << Viewer3D<>::updateDisplay;
            app.exec();
}

int main(int argc, char** argv) {
        testCone(argc, argv);
        return 0;
}
