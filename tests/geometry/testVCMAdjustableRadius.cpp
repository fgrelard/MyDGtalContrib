#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/MedialAxis.h"
#include "geometry/VCMAdjustableRadius.h"
#include "ReverseHatPointFunction.h"

using namespace DGtal;
using namespace std;



void testVCMAdaptable() {
    typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
    typedef Z3i::Object26_6 Graph;
    typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
    typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;

    Image image = GenericReader<Image>::import("/home/fgrelard/test_img/bronche2.vol");
    Z3i::DigitalSet setVolume(image.domain());
    SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

    VCMAdjustableRadius<Z3i::Space, Z3i::L2Metric> vcm(10, ceil( 10 ), Z3i::l2Metric, true );
    vcm.init(setVolume.begin(), setVolume.end());
    DGtal::functors::ReverseHatPointFunction<Z3i::Point, double> chi(1.0, 5.0);
    Ball<Z3i::Point> ball(*setVolume.begin(), 5);
    Z3i::DigitalSet points = ball.pointSet();
    auto mat = vcm.measure(points, chi, *setVolume.begin());
    trace.info() << mat << endl;
    // const Color CURVE3D_COLOR( 100, 100, 140, 128 );
    // for (const Z3i::Point& p : medialAxis.pointSet()) {
    //         viewer << CustomColors3D(Color::Green, Color::Green) << p;
    // }
    // //viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << setVolume;
    // viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    //testBezier(viewer);
    //testBezierDeCasteljau(viewer);
//    testLinking();
//    testDijkstra(viewer);
    testVCMAdaptable();
    app.exec();

    return 0;
}
