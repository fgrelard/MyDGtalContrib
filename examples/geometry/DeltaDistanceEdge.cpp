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

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int argc, char** argv )
{
        using namespace DGtal;
        using namespace DGtal::Z2i;

        typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
        typedef ImageContainerBySTLVector<Domain,float>         FloatImage2D;
        typedef DistanceToMeasureEdge<FloatImage2D>                 Distance;

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


        double min = delta.distance2(*std::min_element(delta.domain().begin(), delta.domain().end(), [&](const Point& f,
                                                                                        const Point& s) {
                                              return delta.distance2(f) < delta.distance2(s);
                                                       }));
        trace.info() << min << endl;
        GradientColorMap<float> cmap_grad( min, m );
        cmap_grad.addColor( Color( 255, 255, 255 ) );
        cmap_grad.addColor( Color( 255, 255, 0 ) );
        cmap_grad.addColor( Color( 255, 0, 0 ) );
        cmap_grad.addColor( Color( 0, 255, 0 ) );
        cmap_grad.addColor( Color( 0,   0, 255 ) );
        cmap_grad.addColor( Color( 0,   0, 0 ) );
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();

        for ( typename Domain::ConstIterator it = delta.domain().begin(),
                      itE = delta.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                Z3i::Point c(p[0], p[1], 0);
                float v = sqrt( delta.distance2( p ) );
                v = std::min( (float)m, std::max( v, 0.0f ) );
                RealVector grad = delta.projection( p );
                Z3i::RealVector grad3D(grad[0], grad[1], 0);
                // board <<  CustomStyle( p.className(),
                //                        new CustomColors( Color::Green, cmap_grad( v ) ) ) <<
                //         Point(p[0] + grad[0], p[1] + grad[1]);
                Color currentColor = cmap_grad( v );
                currentColor.alpha(255);
                viewer << CustomColors3D( Color::Black, currentColor )
                       << c;
                if (grad3D.norm() < 20)
                        viewer.addLine( c, c+grad3D );

        }
        viewer << Viewer3D<>::updateDisplay;
        app.exec();

        return 0;
}
