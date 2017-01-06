#include "shapes/Curve.h"
#include "geometry/CurveProcessor.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "surface/SurfaceUtils.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <QtGui/qapplication.h>

using namespace DGtal;
using namespace std;

void testCurveConnexity(int argc, char** argv) {
        trace.beginBlock("Test curve connexity");
        // QApplication app(argc, argv);
        // Viewer3D<> viewer;
        // viewer.show();

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        Z3i::DigitalSet newCurve = curveProcessor.ensureConnexity();
        trace.info() << newCurve.size() << " " << setVolume.size() << endl;

        // viewer << Viewer3D<>::updateDisplay;
        // app.exec();
         trace.endBlock();
}


void testCurveEndPoints(int argc, char** argv) {
        trace.beginBlock("Test curve endpoints");
        // QApplication app(argc, argv);
        // Viewer3D<> viewer;
        // viewer.show();

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        Z3i::DigitalSet newCurve = curveProcessor.endPoints();
        for (const Z3i::Point& p : newCurve)
                trace.info() << p << endl;

        // viewer << Viewer3D<>::updateDisplay;
        // app.exec();
         trace.endBlock();
}

void testCurveBranchingPoints(int argc, char** argv) {
        trace.beginBlock("Test curve branching points");
        // QApplication app(argc, argv);
        // Viewer3D<> viewer;
        // viewer.show();

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        Z3i::DigitalSet newCurve = curveProcessor.branchingPoints();
        for (const Z3i::Point& p : newCurve)
                trace.info() << p << endl;

        // viewer << Viewer3D<>::updateDisplay;
        // app.exec();
         trace.endBlock();
}

void testCurveOrdered(int argc, char** argv) {
        trace.beginBlock("Test curve ordered");
        // QApplication app(argc, argv);
        // Viewer3D<> viewer;
        // viewer.show();

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/Bezier/skeletonTube250_2.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        Curve<std::vector<Z3i::Point> > newCurve = curveProcessor.convertToOrderedCurve();
        for (const Z3i::Point& p : newCurve.pointSet())
                trace.info() << p << endl;

        // viewer << Viewer3D<>::updateDisplay;
        // app.exec();
         trace.endBlock();
}


int main(int argc, char** argv) {
        testCurveConnexity(argc, argv);
        testCurveEndPoints(argc, argv);
        testCurveOrdered(argc, argv);
        return 0;
}
