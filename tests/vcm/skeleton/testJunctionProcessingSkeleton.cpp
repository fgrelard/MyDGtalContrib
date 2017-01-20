#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <ctime>
#include <vector>
#include <cstdlib>
#include "vcm/skeleton/post/JunctionProcessingSkeleton.h"
#include "geometry/predicate/AbovePlanePredicate.h"

using namespace DGtal;
using namespace std;

void testConstructJunctionProcessing(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef AbovePlanePredicate<Z3i::Space> Predicate;
        typedef DigitalPlane<Z3i::Space> Plane;

        Image image = GenericReader<Image>::import("/home/florent/test_img/junction_simple4PI045.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Image curve = GenericReader<Image>::import("/home/florent/test_img/skeleton20Junction45.vol");
        Z3i::DigitalSet setCurve(curve.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setCurve, curve, 0, 255);
        std::vector<Plane> f = {        Plane({0, 2, 25}, {0, -0.997131, -0.0756897}),
                                        Plane( {0, 37, 25}, {-0, 0.99988, 0.015476}),
                                        Plane( {0, -29, 1}, {0, -0.807488, -0.589884}),
                                        Plane( {0, -4, 21}, {-0, 0.786376, 0.617748}),
                                        Plane( {0, -30, 47}, {0, -0.709414, 0.704792}),
                                        Plane( {0, -28, 49}, {0, -0.737378, 0.675481}),
                                        Plane( {0, -3, 28}, {-0, 0.777478, -0.62891})
        };




        Z3i::DigitalSet a3shell(setCurve.domain());
        a3shell.insert({0, 0, 25});
        a3shell.insert({0, -3, 24});
        a3shell.insert({0, 1, 25});
        a3shell.insert({0, -2, 26});

        for (const Plane& plane : f) {
                viewer << CustomColors3D(Color::Blue, Color::Blue);
                viewer.addLine(plane.getCenter(), plane.getCenter() + plane.getPlaneEquation().normal()*10);
        }


         JunctionProcessingSkeleton<Z3i::DigitalSet, Predicate> junctionProc(setCurve, a3shell, setVolume, f);
        Z3i::DigitalSet link = junctionProc.postProcess();
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << link;

//        Z3i::DigitalSet junctionProcessed = junctionProc.postProcess();
        viewer << Viewer3D<>::updateDisplay;
}


int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testConstructJunctionProcessing(viewer);
        app.exec();
        return 0;
}
