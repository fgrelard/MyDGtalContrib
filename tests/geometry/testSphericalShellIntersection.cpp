#include "geometry/CurveDecomposition.h"
#include "geometry/CurveProcessor.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <QtGui/qapplication.h>
#include "geometry/SphericalShellIntersection.h"

using namespace DGtal;
using namespace std;


void testShell(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/deformed_cylinder.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        Z3i::Point center(0,0,25);
        SphericalShellIntersection<Z3i::DigitalSet> ssi(setVolume, center, 11, 13);
      Z3i::DigitalSet shell = ssi.shell();
      viewer << CustomColors3D(Color::Red, Color::Red) << center;
      viewer << CustomColors3D(Color::Red, Color::Red) << shell;
      viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
}

void testDegree() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        Z3i::Point center(0,0,25);
        SphericalShellIntersection<Z3i::DigitalSet> ssi(setVolume, center, 11, 13);
        unsigned int degree = ssi.degree();
        trace.info() << "Test " << ((degree == 2) ? "passed" : "failed") << endl;
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    testShell(viewer);
    testDegree();
    app.exec();
}
