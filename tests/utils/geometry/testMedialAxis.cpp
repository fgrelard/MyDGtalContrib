#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "geometry/path/LinkPointAlgorithm.h"
#include "geometry/path/BezierCasteljauLinkAlgorithm.h"
#include "geometry/path/BezierLinkAlgorithm.h"
#include "geometry/path/BresenhamAlgorithm.h"
#include "geometry/path/DijkstraAlgorithm.h"
#include "geometry/path/AStarAlgorithm.h"
#include "shapes/Curve.h"
#include "shapes/Border.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/MedialAxis.h"

using namespace DGtal;
using namespace std;



void testMedialAxis(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef Z3i::Object26_6 Graph;
        typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
        typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;

        Image image = GenericReader<Image>::import("/home/florent/test_img/bronche2.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        MedialAxis<Z3i::DigitalSet> medialAxis(setVolume);

        const Color CURVE3D_COLOR( 100, 100, 140, 128 );
        for (const Z3i::Point& p : medialAxis.pointSet()) {
                viewer << CustomColors3D(Color::Green, Color::Green) << p;
        }
        //viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    //testBezier(viewer);
    //testBezierDeCasteljau(viewer);
//    testLinking();
//    testDijkstra(viewer);
    testMedialAxis(viewer);
    app.exec();

    return 0;
}
