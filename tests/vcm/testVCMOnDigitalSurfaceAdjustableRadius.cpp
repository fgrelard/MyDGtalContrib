#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/MedialAxis.h"
#include "vcm/VCMOnDigitalSurfaceAdjustableRadius.h"
#include "ReverseHatPointFunction.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/topology/SurfelAdjacency.h"

using namespace DGtal;
using namespace std;



void testVCMAdaptable() {
    typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
    typedef Z3i::Object26_6 Graph;
    typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
    typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;
    typedef DGtal::ExactPredicateLpSeparableMetric<DGtal::Z3i::Space, 2> Metric;
    typedef DGtal::functors::BallConstantPointFunction<DGtal::Z3i::Point, double> KernelFunction;

    typedef DGtal::Z3i::KSpace KSpace;
    typedef DGtal::ImplicitDigitalSurface< KSpace, Z3i::DigitalSet > DigitalSurfaceContainer;

    typedef KSpace::Surfel Surfel;
    typedef VCMOnDigitalSurfaceAdjustableRadius< DigitalSurfaceContainer, Metric, KernelFunction> VCMOnSurface;

    Image image = GenericReader<Image>::import("/home/fgrelard/test_img/cylinder.vol");
    Z3i::DigitalSet setVolume(image.domain());
    SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

    double R = 10, r = 10, delta = 2;
    Z3i::Domain myDomain = setVolume.domain();
    Metric l2;
    KSpace ks;
    ks.init( myDomain.lowerBound(),
             myDomain.upperBound(), true );
    SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
    Surfel bel = Surfaces<KSpace>::findABel( ks, setVolume, 1000000 );
    DigitalSurfaceContainer* container =
            new DigitalSurfaceContainer( ks, setVolume, surfAdj, bel, false  );
    DigitalSurface<DigitalSurfaceContainer>* mySurface = new DigitalSurface< DigitalSurfaceContainer >( container ); //acquired

    //! [DVCM3D-instantiation]
    Surfel2PointEmbedding embType = Pointels; // Could be Pointels|InnerSpel|OuterSpel;
    KernelFunction chiSurface( 1.0, r );             // hat function with support of radius r

    MedialAxis<DGtal::Z3i::DigitalSet> medialAxis(setVolume);
    VCMOnSurface myVCMSurface( *mySurface, embType, R, r,
                                  chiSurface, medialAxis, delta, l2, true);

    auto matrix = myVCMSurface.mapPoint2ChiVCM().begin()->second.values;
    trace.info() << matrix << endl;
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
    testVCMAdaptable();
    app.exec();

    return 0;
}
