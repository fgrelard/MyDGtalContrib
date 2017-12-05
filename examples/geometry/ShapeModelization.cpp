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


Z3i::DigitalSet translate(const Z3i::DigitalSet& toTranslate,
                          const Z3i::RealVector& v) {

    Z3i::DigitalSet translated(toTranslate.domain());
    for (const Z3i::Point&  p :toTranslate) {
        translated.insert(p+v);
    }
    return translated;
}

int main( int argc, char** argv )
{
    typedef DGtal::DistanceTransformation<Z3i::Space, Z3i::DigitalSet, Z3i::L2Metric> DTL2;
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

    Ball<Z3i::Point> ball(Z3i::Point(0,0,0), 5);
    Z3i::DigitalSet ballPoints = ball.pointSet();
    Z3i::DigitalSet bt1 = translate(ballPoints, Z3i::RealVector(6,6,6));
    Z3i::DigitalSet bt2 = translate(ballPoints, Z3i::RealVector(6,18,6));
    Z3i::DigitalSet bt3 = translate(ballPoints, Z3i::RealVector(6,30,6));

    Z3i::DigitalSet cylinder = modeller.drawCylinder(30, 5);
    Z3i::DigitalSet ct1 = translate(cylinder, Z3i::RealVector(30,6,6));
    Z3i::DigitalSet ct2 = translate(cylinder, Z3i::RealVector(30,16,6));

    Z3i::DigitalSet ct3 = translate(cylinder, Z3i::RealVector(70,6, 6));
    Z3i::DigitalSet ct4 = translate(cylinder, Z3i::RealVector(70,26,6));
    Z3i::DigitalSet surface = modeller.drawSurface(30);
    Z3i::DigitalSet st = translate(surface, Z3i::RealVector(50,50,50));

    bt1.insert(bt2.begin(), bt2.end());
    bt1.insert(bt3.begin(), bt3.end());

    bt1.insert(ct1.begin(), ct1.end());
    bt1.insert(ct2.begin(), ct2.end());

    bt1.insert(ct3.begin(), ct3.end());
    bt1.insert(ct4.begin(), ct4.end());

    bt1.insert(st.begin(), st.end());

    Z3i::Domain domain(Z3i::Point(-100,-100,-100)-Z3i::Point::diagonal(2), Z3i::Point(100, 300, 300)+Z3i::Point::diagonal(2));


    ct1.insert(ct2.begin(), ct2.end());
    ct1.insert(ct3.begin(), ct3.end());
    ct1.insert(ct4.begin(), ct4.end());

    Point lower, upper;
    Z3i::DigitalSet airway(Z3i::Domain(Z3i::Point(0,0,0), Z3i::Point(100,100,100)));

    //modeller.createDeformedSyntheticAirwayTree(airway, 4, 50, 0, 0, {0,0,0});

    airway = modeller.createHelixCurve(20,10,5);
    airway.computeBoundingBox(lower, upper);


    Z3i::Domain domainObjects(lower - Point::diagonal(1), upper + Point::diagonal(1));
    Z3i::L2Metric l2Metric;
    DTL2 dt(domainObjects, airway, l2Metric);

    Image3D anImage3D(domainObjects);
    for (const Z3i::Point& p : domainObjects) {
        anImage3D.setValue(p, 120);
    }
    for (auto it = domainObjects.begin(), ite = domainObjects.end();
         it != ite; ++it) {
        if (airway.find(*it) != airway.end())
            anImage3D.setValue(*it, 255);
            //anImage3D.setValue(*it, 140 + dt(*it) * 10);
    }
//	anImage3D = ImageFromSet<Image3D>::create(set, 1);

    VolWriter<Image3D>::exportVol(outfile, anImage3D);
    const Color CURVE3D_COLOR( 100, 100, 140, 128 );


    trace.info() << "saved" << endl;
}
///////////////////////////////////////////////////////////////////////////////
