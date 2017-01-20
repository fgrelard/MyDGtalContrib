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
#include "geometry/CurveDecomposition.h"
#include "vcm/CuttingPlaneEstimator.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "vcm/RecenterSkeleton.h"

using namespace DGtal;
using namespace std;

void testRecentering(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef functors::BallConstantPointFunction<Z3i::Point, double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> PlaneEstimator;
        typedef DistanceTransformation<Z3i::Space, Z3i::DigitalSet, Z3i::L2Metric> DTL2;


        Image image = GenericReader<Image>::import("/home/florent/test_img/junction_simple4PI045.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Image curve = GenericReader<Image>::import("/home/florent/test_img/Bezier/skeletonJunctionPI045Connexity.vol");
        Z3i::DigitalSet setCurve(curve.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setCurve, curve, 0, 255);

        RecenterSkeleton<Z3i::DigitalSet> recentering(setCurve, setVolume);
        Z3i::DigitalSet skeletonRecentered = recentering.recenter();
        viewer << CustomColors3D(Color::Blue, Color::Blue) << skeletonRecentered;
        viewer << CustomColors3D(Color(210,210,210,20), Color(210,210,210,20)) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
        srand(time(NULL));
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testRecentering(viewer);
        app.exec();
        return 0;
}
