#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/MedialAxis.h"
#include "shapes/Border.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "shapes/DigitalPlane.h"
#include "geometry/DigitalPlaneProcessor.h"
#include "ShapeDescriptor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


Z3i::DigitalSet
markPointsBetweenPlanes(const Z3i::DigitalSet& setVolume,
                        const DigitalPlane<Z3i::Space> &currentPlane,
                        const DigitalPlane<Z3i::Space> &previousPlane,
                        double distanceMax) {

    typedef DigitalPlane<Z3i::Space> Plane;
    Z3i::L2Metric l2Metric;
    Z3i::DigitalSet difference(setVolume.domain());
    Z3i::RealVector dirCurrent = currentPlane.getPlaneEquation().normal();
    Z3i::RealVector dirPrevious = previousPlane.getPlaneEquation().normal();
    Z3i::Point pCurrent = currentPlane.getCenter();
    Z3i::Point pPrevious = previousPlane.getCenter();

    if (dirPrevious.dot(dirCurrent) >= 0) {
        dirCurrent = -dirCurrent;
    }

    Plane currentPlane2(pCurrent, dirCurrent, currentPlane.getConnexity());

    for (const Z3i::Point &p : setVolume) {
        if (l2Metric(p, pCurrent) > distanceMax) continue;
        if (currentPlane2.isPointAbove(p) &&
            previousPlane.isPointAbove(p))
            difference.insert(p);
    }
    return difference;
}


int main( int  argc, char**  argv )
{
    typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
    typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
    typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> OrthoPlaneEstimator;

    typedef DGtal::ExactPredicateLpSeparableMetric<Z3i::Space, 2> L2Metric;
    typedef DGtal::DistanceTransformation<Z3i::Space, Z3i::DigitalSet, L2Metric> DTL2;
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
        ("thresholdMin,m", po::value<int>()->default_value(1), "minimum threshold for binarization")
        ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
        ("step,z", po::value<int>()->default_value(10), "z step")
        ;

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    } catch(const std::exception& ex){
        parseOK=false;
        trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
    po::notify(vm);
    if( !parseOK || vm.count("help")||argc<=1)
    {
        std::cout << "Usage: " << argv[0] << " [input]\n"
                  << "Display volume file as a voxel set by using QGLviewer"<< endl
                  << general_opt << "\n";
        return 0;
    }
    if(!vm.count("input"))
    {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }

    string inputFilename = vm["input"].as<std::string>();
    int thresholdMin = vm["thresholdMin"].as<int>();
    int thresholdMax = vm["thresholdMax"].as<int>();
    int z = vm["step"].as<int>();
    Image volume = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domainVolume = volume.domain();
    Z3i::DigitalSet setVolume(domainVolume);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                  thresholdMin-1, thresholdMax);
    QApplication application(argc,argv);
    Viewer3D<> viewer;
    viewer.show();

    MedialAxis<Z3i::DigitalSet> maComputer(setVolume);
    Z3i::DigitalSet ma = maComputer.pointSet();
    Border<Z3i::DigitalSet> borderComputer(setVolume);
    Z3i::DigitalSet border = borderComputer.pointSet();
    Z3i::Point p = *std::next(ma.begin(), z);
    DGtal::trace.info() << p << std::endl;
    Z3i::DigitalSet geodesicLoop(border.domain());
    Z3i::DigitalSet proj(border.domain());
    SetProcessor<Z3i::DigitalSet> setProc(border);
    Z3i::Point pSurface = setProc.closestPointAt(p);
    proj.insert(pSurface);
    Z3i::Object26_6 obj(Z3i::dt26_6, setVolume);
    std::vector<Z3i::Point> neighbors;

    std::back_insert_iterator<std::vector<Z3i::Point> > inserter(neighbors);
    obj.writeNeighbors(inserter, p);
    for (const Z3i::Point& n : neighbors) {
        Z3i::Point pSurface = setProc.closestPointAt(n);
        proj.insert(pSurface);
    }

    for (const Z3i::Point& f : proj) {
        for (const Z3i::Point& s : proj) {
            if (f == s) continue;
            AStarAlgorithm<Z3i::Point, Z3i::DigitalSet> aStar(f, s, border);
            std::vector<Z3i::Point> geodesic = aStar.linkPoints();
            geodesicLoop.insert(geodesic.begin(), geodesic.end());
        }
    }

    std::vector<Z3i::Point> loops(geodesicLoop.begin(), geodesicLoop.end());
        ShapeDescriptor<Z3i::DigitalSet> stats(geodesicLoop);
    Z3i::RealVector normal = stats.computeNormalFromLinearRegression();
    DigitalPlane<Z3i::Space> plane(p, normal);
    DigitalPlaneProcessor<Z3i::Space> digPlaneProc(plane);
    std::vector<Z3i::RealVector> points = digPlaneProc.planeToQuadrangle();
    double f = 20.0;
    viewer.setLineColor(Color::Blue);
    viewer.setFillColor(Color::Blue);
    viewer.setFillTransparency(150);
    viewer.addQuad(p + points[0] * f,
                   p + points[1] * f,
                   p + points[2] * f,
                   p + points[3] * f);

    //viewer << CustomColors3D(Color::Yellow, Color::Yellow) << proj;

    for (const Z3i::Point& g : geodesicLoop) {
        viewer << CustomColors3D(Color::Red, Color::Red) << g;
    }




    viewer << CustomColors3D(Color(220,220,220, 50), Color(220,220,220,50)) << setVolume;
    viewer << Viewer3D<>::updateDisplay;
    application.exec();
    return 0;


}
