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
    Z3i::DigitalSet deformedcylinder = modeller.drawDeformedCylinder(50, 10);
    Z3i::DigitalSet translated(cylinder.domain());
    for (const auto& p : deformedcylinder)
        translated.insert(p + Z3i::Point(2, 0, 0));
    cylinder.insert(translated.begin(), translated.end());
    Z3i::Domain domain(Z3i::Point(-100,-100,-100)-Z3i::Point::diagonal(2), Z3i::Point(100, 300, 300)+Z3i::Point::diagonal(2));


    Point lower, upper;

    cylinder.computeBoundingBox(lower, upper);
    Z3i::Domain domainObjects(lower - Point::diagonal(1), upper + Point::diagonal(1));
    Image3D anImage3D(domainObjects);
    for (auto it = domainObjects.begin(), ite = domainObjects.end();
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
