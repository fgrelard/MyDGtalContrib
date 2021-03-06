#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "geometry/DistanceToMeasureEdge.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "vcm/DeltaVCM.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

namespace DGtal {
    typedef SimpleMatrix<double, 2, 2> MatrixDouble;

    bool operator!=(const MatrixDouble &m1, const MatrixDouble &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<float, 2, 2> MatrixFloat;

    bool operator!=(const MatrixFloat &m1, const MatrixFloat &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<double, 3, 3> MatrixDouble3D;

    bool operator!=(const MatrixDouble3D &m1, const MatrixDouble3D &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<float, 3, 3> MatrixFloat3D;

    bool operator!=(const MatrixFloat3D &m1, const MatrixFloat3D &m2) { return !(m1 == m2); }

    namespace functors {
        bool operator==(Identity f1, Identity f2) { return true; }
    }
}

int main(int argc, char **argv) {
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
    typedef ImageContainerBySTLVector<Domain,double>         FloatImage2D;
    typedef DistanceToMeasureEdge<FloatImage2D>                 Distance;
    typedef DeltaVCM<Distance> DVCM;
    typedef functors::BallConstantPointFunction<Z2i::Point, double> KernelFunction;
    typedef EigenDecomposition<2,double> LinearAlgebraTool;
    typedef PointVector<2,double> RealVector2f;
    typedef typename DVCM::Matrix Matrix;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
            ("mass,a", po::value<double>()->default_value(1), "Mass to integrate for distance to measure")
            ("rmax,r", po::value<double>()->default_value(10), "Max radius for delta distance")
            ("radiusInside,R", po::value<double>()->default_value(10), "radius of the ball inside voronoi cell")
            ("radiusNeighbour,n", po::value<double>()->default_value(10), "radius of the ball for the neighbourhood")
            ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
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
    int thresholdMax = vm["thresholdMax"].as<int>();
    double mass = vm["mass"].as<double>();
    double rmax = vm["rmax"].as<double>();
    double R = vm["radiusInside"].as<double>();
    double r = vm["radiusNeighbour"].as<double>();

    GrayLevelImage2D img  = GenericReader<GrayLevelImage2D>::import( inputFilename );
    auto domain = img.domain();
    FloatImage2D     fimg( img.domain() );
    FloatImage2D::Iterator outIt = fimg.begin();
    for ( GrayLevelImage2D::ConstIterator it = img.begin(), itE = img.end();
          it != itE; ++it )
    {
        float v = ((float)*it) * 1.0 / thresholdMax;
        *outIt++ = v;
    }
    trace.beginBlock( "Computing delta-distance." );
    Distance     delta( mass, fimg, rmax );
    trace.endBlock();

    KernelFunction chi( 1.0, r );
    DVCM dvcm(delta, R, r);

    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();

    Matrix vcm_r, evec, nil;
    RealVector2f eval;

    std::map<Point, double> pToTube;
    for ( typename Domain::ConstIterator it = delta.domain().begin(),
                  itE = delta.domain().end(); it != itE; ++it ) {

        Point p = *it;
        vcm_r = dvcm.measure(chi, p);
        if (vcm_r == nil) continue;
        LinearAlgebraTool::getEigenDecomposition(vcm_r, evec, eval);
        double l0 = eval[0];
        double l1 = eval[1];
        double tub = l1 / (l0 + l1);
        pToTube[p] = tub;
    }

    double min = std::min_element(pToTube.begin(), pToTube.end(), [&](const std::pair<Point, double>& p1,
                                                                      const std::pair<Point, double>& p2) {
        return p1.second < p2.second;
    })->second;
    double max = std::max_element(pToTube.begin(), pToTube.end(), [&](const std::pair<Point, double>& p1,
                                                                      const std::pair<Point, double>& p2) {
        return p1.second < p2.second;
    })->second;

    GradientColorMap<float> cmap_grad( min, max );
    cmap_grad.addColor( Color( 255, 255, 255 ) );
    cmap_grad.addColor( Color( 255, 255, 0 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 0, 255, 0 ) );
    cmap_grad.addColor( Color( 0,   0, 255 ) );
    cmap_grad.addColor( Color( 0,   0, 0 ) );

    int index = 0;
    for (const auto& pair : pToTube) {
        Point p = pair.first;
        Z3i::Point p3d(p[0], p[1], 0);
        double val = pair.second;
        Color currentColor = cmap_grad(val);
        viewer << CustomColors3D(currentColor, currentColor) << p3d;
        index++;
    }

    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
