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
#include "geometry/DistanceToMeasure.h"
#include "shapes/Ball.h"
#include "geometry/SetProcessor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

double signedAngle(const Z2i::RealVector& v) {
        double angle = atan2(v[1], v[0]);
        return (angle > 0 ? angle : (2*M_PI + angle));
}

double signedAngleTwoVectors(const Z2i::RealVector& v1,
                             const Z2i::RealVector& v2) {
        if (v1 == v2) return 2 * M_PI;
        double angle = signedAngle(v1) - signedAngle(v2);
        if (angle < 0) {
                return 2 * M_PI;
        }
        return angle;
}

double twoAngleClusters(const std::vector<Z2i::RealVector>& dirs) {
        if (dirs.size() == 1) return 0;
        double maxAngle = 0;
        for (const Z2i::RealVector& f : dirs) {
                Z2i::RealVector other = *min_element(dirs.begin(),
                                                     dirs.end(),
                                                     [&](const Z2i::RealVector& lhs,
                                                         const Z2i::RealVector& rhs) {
                                                             return (signedAngleTwoVectors(f, lhs) < signedAngleTwoVectors(f, rhs));
                        });

                double angle = signedAngleTwoVectors(f, other);
                if (angle != 2* M_PI && angle > maxAngle)
                        maxAngle = angle;
        }
        return maxAngle;
}

int main( int argc, char** argv )
{
        using namespace DGtal;
        using namespace DGtal::Z2i;

        typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
        typedef ImageContainerBySTLVector<Domain,float>         FloatImage2D;
        typedef DistanceToMeasure<FloatImage2D>                 Distance;

        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
                ("mass,a", po::value<double>()->default_value(1), "Mass to integrate for distance to measure")
                ("rmax,r", po::value<double>()->default_value(10), "Max radius for delta distance")
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

        float m = 0.0f;
        for ( typename Domain::ConstIterator it = delta.domain().begin(),
                      itE = delta.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                float v = sqrt( delta.distance2( p ) );
                m = std::max( v, m );
        }

        //  GradientColorMap<float> cmap_grad( 0, 1.82 );
        GradientColorMap<float> cmap_grad( 0, 2*M_PI );

        cmap_grad.addColor( Color( 255, 255, 255 ) );
        cmap_grad.addColor( Color( 255, 255, 0 ) );
        cmap_grad.addColor( Color( 255, 0, 0 ) );
        cmap_grad.addColor( Color( 0, 255, 0 ) );
        cmap_grad.addColor( Color( 0,   0, 255 ) );
        cmap_grad.addColor( Color( 0,   0, 0 ) );
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();



        vector<double> ve;
        for ( typename Domain::ConstIterator it = delta.domain().begin(),
                      itE = delta.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                Z3i::Point c(p[0], p[1], 0);
                float v = sqrt( delta.distance2( p ) );
                v = std::min( (float)m, std::max( v, 0.0f ) );
                int radius = 4;
                Ball<Point> ball(p, radius);
                FloatImage2D aSet = ball.intersection( fimg );
                std::vector<RealVector> dirs;
                DigitalSet projections(img.domain());
                for (const Point& b : aSet.domain()) {
                        if (aSet(b) == 0) continue;
                        RealVector grad = delta.projection(b);
                        if (ball.contains(b+grad) &&
                            grad != RealVector::zero) {
                                grad = grad.getNormalized();
                                 dirs.push_back(grad);
                         }
                }

                double angle = twoAngleClusters(dirs);
                RealVector grad = delta.projection(p);
                Z3i::RealVector grad3D(grad[0], grad[1], 0);
                // if (distanceProj < v) continue;
                Color currentColor = cmap_grad(angle);
//                            Color currentColor = cmap_grad( angle );
                currentColor.alpha(255);
                viewer << CustomColors3D( Color::Black, currentColor )
                       << c;
                viewer.addLine( c, c+grad3D*5 );



        }

        viewer << Viewer3D<>::updateDisplay;
        app.exec();

        return 0;
}
