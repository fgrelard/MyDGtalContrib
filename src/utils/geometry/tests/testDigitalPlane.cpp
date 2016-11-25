#include "../DigitalPlane.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "surface/SurfaceUtils.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <QtGui/qapplication.h>

using namespace std;
using namespace DGtal;

void testDigitalPlane(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();

        Z3i::RealVector normal(0.0486464, -0.750268, 0.659342);
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/Pruning/Thinvox/bronche2_512_512_538.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
         Z3i::Point center(159, 228, 170);

        DigitalPlane<Z3i::Space> plane2(center, normal,  6);
        Z3i::DigitalSet planeSet2 = plane2.intersectionWithSetOneCC(setVolume);
        viewer << CustomColors3D(Color::Green, Color::Green) << center;
        viewer << CustomColors3D(Color::Red, Color::Red) <<planeSet2;
        viewer << CustomColors3D(Color(128,128,128,20), Color(128,128,128,20)) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
        app.exec();

}

int main(int argc, char** argv) {
        testDigitalPlane(argc, argv);
        return 0;
}
