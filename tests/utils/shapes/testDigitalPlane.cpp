#include "shapes/DigitalPlane.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <QtGui/qapplication.h>

using namespace std;
using namespace DGtal;

void testEquation() {
        trace.beginBlock("Test equation");
        Z3i::RealVector normal(0.0486464, -0.750268, 0.659342);
        Z3i::Point center(159, 228, 170);

        DigitalPlane<Z3i::Space> plane2(center, normal,  6);
        int expectedD = (normal[0] * center[0] + normal[1] * center[1] + normal[2] * center[2]) * 1000;
        int expectedOmega = (std::abs(normal[0]) + std::abs(normal[1]) + std::abs(normal[2])) * 1000;
        if (plane2.getCenter() == center &&
            plane2.getPlaneEquation().normal() == normal.getNormalized() &&
            expectedD == (int)(plane2.getPlaneEquation().mu()*1000) &&
            expectedOmega == (int)(plane2.getPlaneEquation().nu()*1000))
                trace.info() << "Passed" << endl;
        else
                trace.info() << "Failed" << endl;
        trace.endBlock();
}

void testContains() {
        trace.beginBlock("Test contains");
        Z3i::RealVector normal(0.0486464, -0.750268, 0.659342);
        Z3i::Point center(159, 228, 170);

        DigitalPlane<Z3i::Space> plane2(center, normal,  6);
        Z3i::Point other(160,228,170);
        if (plane2.contains(center) && plane2.contains(other))
                trace.info() << "Passed" << endl;
        else
                trace.info() << "Failed" << endl;
        trace.endBlock();
}

void testPointAbove() {
        trace.beginBlock("Test point above");
        Z3i::RealVector normal(0.0486464, -0.750268, 0.659342);
        Z3i::Point center(159, 228, 170);

        DigitalPlane<Z3i::Space> plane2(center, normal,  6);
        Z3i::Point other(158,227,169);
        if (plane2.isPointAbove(center) && plane2.isPointAbove(other))
                trace.info() << "Passed" << endl;
        else
                trace.info() << "Failed" << endl;
        trace.endBlock();

}

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
        testEquation();
        testContains();
        testPointAbove();
        testDigitalPlane(argc, argv);
        return 0;
}
