#include "geometry/CurveProcessor.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <QtGui/qapplication.h>
#include <iterator>

using namespace DGtal;
using namespace std;

void testCurveConnectivity() {
        trace.beginBlock("Test curve connectivity");

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        Z3i::DigitalSet newCurve = curveProcessor.ensureConnectivity();
        trace.info() << newCurve.size() << " " << setVolume.size() << endl;

         trace.endBlock();
}


void testCurveEndPoints() {
        trace.beginBlock("Test curve endpoints");

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        Z3i::DigitalSet newCurve = curveProcessor.endPoints();

        for (const Z3i::Point& p : newCurve)
                trace.info() << p << endl;
         trace.endBlock();
}

void testCurveBranchingPoints() {
        trace.beginBlock("Test curve branching points");


        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);
        Z3i::DigitalSet newCurve = curveProcessor.branchingPoints();
        for (const Z3i::Point& p : newCurve)
                trace.info() << p << endl;

         trace.endBlock();
}

void testCurveOrdered() {
        trace.beginBlock("Test curve ordered");

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/Bezier/skeletonTube250_2.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);

        std::vector<Z3i::Point> newCurve = curveProcessor.convertToOrderedCurve();
        for (const Z3i::Point& p : newCurve)
                trace.info() << p << endl;

         trace.endBlock();
}

void testCurveEnsureOneCC(Viewer3D<>& viewer) {
           typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/airway.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Image curve = GenericReader<Image>::import("/home/florent/test_img/skeletonVCMAirway.vol");
        Z3i::DigitalSet setCurve(curve.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setCurve, curve, 0, 255);
        CurveProcessor<Z3i::DigitalSet> curveProcessor(setCurve);
        double min = std::numeric_limits<double>::min();
        double max = std::numeric_limits<double>::max();
        Z3i::DigitalSet curveCC = curveProcessor.fillHoles(setVolume);
        Z3i::Object26_6 obj(Z3i::dt26_6, curveCC);
        vector<Z3i::Object26_6> cc;
        std::back_insert_iterator< vector<Z3i::Object26_6> > inserter(cc);
        unsigned int nbCC = obj.writeComponents(inserter);
        trace.info() << "Test " << ((nbCC == 1) ? "passed" : " failed") << endl;
        viewer << CustomColors3D(Color::Red, Color::Red) << setCurve;
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << curveCC;
}

void testCurveEnsureOneCCDistance(Viewer3D<>& viewer) {
           typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);


        Z3i::DigitalSet setCurve(image.domain());
        setCurve.insert(Z3i::Point(0,0,5));
        setCurve.insert(Z3i::Point(0,1,5));
        setCurve.insert(Z3i::Point(0,0,7));
        setCurve.insert(Z3i::Point(0,1,9));
        setCurve.insert(Z3i::Point(0,3,11));
        setCurve.insert(Z3i::Point(2,5,13));
        setCurve.insert(Z3i::Point(2,5,17));

        CurveProcessor<Z3i::DigitalSet> curveProcessor(setCurve);
        Z3i::DigitalSet curveCC = curveProcessor.fillHoles( sqrt(3), 2 * sqrt(3));
        Z3i::Object26_6 obj(Z3i::dt26_6, curveCC);
        vector<Z3i::Object26_6> cc;
        std::back_insert_iterator< vector<Z3i::Object26_6> > inserter(cc);
        unsigned int nbCC = obj.writeComponents(inserter);
        viewer << CustomColors3D(Color::Red, Color::Red) << setCurve;
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << curveCC;
}


int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testCurveConnectivity();
        testCurveEndPoints();
        testCurveOrdered();
        //testCurveEnsureOneCC(viewer);
        testCurveEnsureOneCCDistance(viewer);
        viewer << Viewer3D<>::updateDisplay;
        app.exec();
        return 0;
}
