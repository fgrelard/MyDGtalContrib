/** File allowing to create a curve in DGtal according to an equation **/
#include <iostream>
#include <fstream>
#include <time.h>
#include "DGtal/base/Common.h"
#include "DGtal/io/Display3DFactory.h"
// Shape construction
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
// Drawing
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/io/writers/VolWriter.h"

#include "geometry/Distance.h"
#include "shapes/Ball.h"
#include "Modeller.h"
///////////////////////////////////////////////////////////////////////////////
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/io/boards/Board2D.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template <typename Point>
class PointDeletable : public Point {
public:
    typedef typename Point::Component Scalar;
public:
    PointDeletable() : Point() {}
    PointDeletable(const Scalar& x, const Scalar& y, const Scalar& z) : Point(x,y,z) {}
    bool myMarkedToDelete = false;
};


typedef Z3i::DigitalSet DigitalSet;
typedef Z3i::Space Space;
typedef Z3i::KSpace KSpace;
typedef PointDeletable<Z3i::Point> Point;
typedef vector<Point>::const_iterator PointIterator;


int main( int argc, char** argv )
{
    srand(time(NULL));
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("output,o", po::value<std::string>(), "output file")
            ("increment,i", po::value<float>()->default_value(0.01), "increment")
            ("noise,n", po::value<float>()->default_value(0.0f), "noise proba (0-1)")
            ;

    bool parseOK = true;
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, (const char *const *) argv, general_opt), vm);
    } catch (const std::exception &ex) {
        parseOK = false;
        trace.info() << "Error checking program options: " << ex.what() << endl;
    }
    po::notify(vm);
    if (!parseOK || vm.count("help") || argc <= 1) {
        std::cout << "Usage: " << argv[0] << " input output [options]\n"
                  << "Compute the Frangi Vesselness Filter" << endl
                  << general_opt << "\n";
        return 0;
    }
    if (!vm.count("output")) {
        trace.error() << " Missing filename(s) " << endl;
        return 0;
    }
    string outfile = vm["output"].as<std::string>();
    float increment = vm["increment"].as<float>();
    float noise = vm["noise"].as<float>();

    typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3D;
    typedef ImageFromSet<Image3D> SetConverter;

    Modeller<Z3i::DigitalSet> modeller;

    int range = 200;
    int pitch =  20;
    int radius = 10;

    Z3i::DigitalSet cylinder = modeller.drawCylinder(50, 10);
    DigitalPlane<Z3i::Space> plane(Z3i::Point(0,75,0), Z3i::RealVector(0,1,0));
    Z3i::Domain domainPlane(Z3i::Point(-5,74,0), Z3i::Point(20,76,50));
    Z3i::DigitalSet aSet(domainPlane);
    aSet.insert(domainPlane.begin(), domainPlane.end());
    Z3i::DigitalSet planeSet = plane.intersectionWithSet(aSet);

    Ball<Z3i::Point> ball(Z3i::Point(0,40,35), 12);
    Z3i::DigitalSet ballSet = ball.pointSet();
    cylinder.insert(ballSet.begin(), ballSet.end());
    cylinder.insert(planeSet.begin(), planeSet.end());
//	drawCone(curve, 32, 15, increment);

//	createContinuousLogarithmicCurve(curve, 50, increment);
//	construct26ConnectedCurve(curve);
//	createVolumeFromCurve(curve, vectorPoints, 7);

//	createRotatedVolumeFromCurve(curve, vectorPoints, 5, 2*M_PI/3);

//	createJunction(curve, vectorPoints, 0.5);
//	Ball<PointVector<3, double>> ball(Point(0,0,0), 10);
    Z3i::Domain domain(Z3i::Point(-100,-100,-100)-Z3i::Point::diagonal(2), Z3i::Point(100, 300, 300)+Z3i::Point::diagonal(2));
    //domain = Z3i::Domain(Z3i::Point(-20,-20,-20), Z3i::Point(20,20,60));
//	createVolumeFromCurve(curve, vectorPoints, 10);
    //createHelixCurve(vectorPoints, range, radius, pitch, increment);
//	drawCircle(vectorPoints, 50.0, 0., 0., 0., increment);
//	createSyntheticAirwayTree(vectorPoints, 4, 40, 0, 0, {10,50,0}, increment);
//	Image3D anImage3D(domain);

//	create2DNaive();
//	vectorPoints = ball.pointsInBall();
//    DigitalSet set(domain);
//    for (auto it = vectorPoints.begin(), itE = vectorPoints.end(); it != itE; ++it) {
//        set.insert(*it);
//    }
//    DigitalSet set2 = addNoise(set, noise);
//	imageFromRangeAndValue(curve.begin(), curve.end(), anImage3D, 150);

    Image3D anImage3D(domain);
    for (auto it = domain.begin(), ite = domain.end();
         it != ite; ++it) {
        if (cylinder.find(*it) != cylinder.end())
            anImage3D.setValue(*it, 255);
    }
//	anImage3D = ImageFromSet<Image3D>::create(set, 1);

    VolWriter<Image3D>::exportVol(outfile, anImage3D);
    const Color CURVE3D_COLOR( 100, 100, 140, 128 );


    trace.info() << "saved" << endl;
}
///////////////////////////////////////////////////////////////////////////////
