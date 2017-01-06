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
#include "geometry/DigitalPlaneProcessor.h"
#include "shapes/DigitalPlane.h"

using namespace DGtal;
using namespace std;


void testPlaneQuadrangle() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef Z3i::Object26_6 Graph;
        typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
        typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;

        Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");

        Z3i::Point center(0,0,5);
        Z3i::RealVector normal(0.72, 0, 0.28);
        DigitalPlane<Z3i::Space> digPlane(center, normal.getNormalized(),  6);
        DigitalPlaneProcessor<Z3i::Space> planeProcessor(digPlane);

        std::vector<Z3i::RealVector> plane = planeProcessor.planeToQuadrangle();
        trace.info() << plane[0] << " " << plane[1] << " " << plane[2] << " " << plane[3] << endl;
}

void testSlice() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image2D;
        typedef Z3i::Object26_6 Graph;
        typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
        typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;

        Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");

        Z3i::Point center(0,0,5);
        Z3i::RealVector normal(0, 0, 1);
        DigitalPlane<Z3i::Space> digPlane(center, normal.getNormalized(),  6);
        DigitalPlaneProcessor<Z3i::Space> planeProcessor(digPlane);

        Image2D slice = planeProcessor.sliceFromPlane(image, 50);
        PGMWriter<Image2D>::exportPGM("slice.pgm", slice);
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    testPlaneQuadrangle();
    testSlice();
    app.exec();

    return 0;
}
