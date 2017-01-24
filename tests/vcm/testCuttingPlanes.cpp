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

using namespace DGtal;
using namespace std;

void testCuttingPlanes(Viewer3D<>& viewer) {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef functors::BallConstantPointFunction<Z3i::Point, double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> PlaneEstimator;
        typedef DistanceTransformation<Z3i::Space, Z3i::DigitalSet, Z3i::L2Metric> DTL2;
        typedef DigitalPlane<Z3i::Space> Plane;


        Image image = GenericReader<Image>::import("/home/florent/test_img/junction_simple4PI045.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Image curve = GenericReader<Image>::import("/home/florent/test_img/Bezier/skeletonJunctionPI045Connexity.vol");
        Z3i::DigitalSet setCurve(curve.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setCurve, curve, 0, 255);

        DTL2 dt(setVolume.domain(), setVolume, Z3i::l2Metric);
        Z3i::DigitalSet branching = CurveProcessor<Z3i::DigitalSet>(setCurve).branchingPoints();

        std::vector<Z3i::DigitalSet> branches;
        auto decompo = CurveDecomposition<Z3i::DigitalSet>(setCurve, branching).branchDecomposition();

        Z3i::Point b = {0,1,25};
        Ball<Z3i::Point> ball(b, 10);
        for (const Z3i::DigitalSet& b : decompo) {
                if (b.size() <= 2) continue;
                Z3i::DigitalSet intersection = ball.intersection(b);
                branches.push_back(intersection);
        }
        KernelFunction chi(1.0, 10);
        PlaneEstimator estimator(setVolume, chi, 20, 10, 6);
        CuttingPlaneEstimator<Plane> cuttingPlaneEstimator(branches, setVolume);
        std::vector<PlaneEstimator::Plane> planes = cuttingPlaneEstimator.cuttingPlanes(estimator, dt);
        for (const PlaneEstimator::Plane& plane : planes) {
                //trace.info() << plane.getCenter() << endl;
                //viewer << CustomColors3D(Color::Yellow, Color::Yellow) << plane.getCenter();
                viewer << CustomColors3D(Color::Red, Color::Red) << plane.intersectionWithSetOneCC(setVolume);
        }
        for (const Z3i::DigitalSet& b : decompo) {
                if (b.size() <= 2) continue;
                int r = rand() % 256, g = rand() % 256, bl = rand() %  256;
                viewer << CustomColors3D(Color(r,g,bl), Color(r,g,bl)) << b;
                branches.push_back(b);
        }
        viewer << CustomColors3D(Color::Blue, Color::Blue) << setCurve;
        viewer << CustomColors3D(Color(210,210,210,20), Color(210,210,210,20)) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
}

int main(int argc, char** argv) {
        srand(time(NULL));
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testCuttingPlanes(viewer);
        app.exec();
        return 0;
}
