#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/PointUtil.h"

using namespace DGtal;
using namespace std;


void testLineTracking() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;

        std::string home = getenv("HOME");
        std::string path = home + "/test_img/cylinder.vol";
        Image image = GenericReader<Image>::import(path);
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        Z3i::Point point(0,0,20);
        Z3i::RealVector dir(0,1,0);
        Z3i::DigitalSet traversedLine = PointUtil::traversedLineTracking(setVolume, point, dir);
        for (const Z3i::Point& p : traversedLine)
                trace.info() << p << " ";
        trace.info() << endl;

}

void testClosestPointTrackingWithNormal(Viewer3D<> & viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image2D;
        typedef Z3i::Object26_6 Graph;
        std::string home = getenv("HOME");
        std::string path = home + "/test_img/cylinder.vol";
        Image image = GenericReader<Image>::import(path);
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Z3i::Point center(0,0,5);
        Z3i::RealVector normal(0, 0, 1);
        Z3i::Point other(3,3,20);
        Z3i::RealVector normalOther(-0.1,-0.3,-0.7);
        std::pair<Z3i::Point, Z3i::Point> pair = PointUtil::twoClosestPointsTrackingWithNormal(setVolume, center, normal, other, normalOther);
        trace.info() << pair.first << " " << pair.second << std::endl;
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << center << other;
        viewer << CustomColors3D(Color::Red, Color::Red) << pair.first << pair.second;
        viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testLineTracking();
        testClosestPointTrackingWithNormal(viewer);
        app.exec();
        return 0;
}
