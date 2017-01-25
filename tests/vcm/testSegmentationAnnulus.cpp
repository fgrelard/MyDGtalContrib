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

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int argc, char** argv )
{
        using namespace DGtal;
        using namespace DGtal::Z2i;

        typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
        typedef ImageContainerBySTLVector<Domain,float>         FloatImage2D;
        typedef DistanceToMeasure<FloatImage2D>                 Distance;
        if ( argc <= 3 ) return 1;

        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
                ("output,o", po::value<std::string>(), "vol file (corresponding volume")
                ("thresholdPruning,t", po::value<double>()->default_value(25), "threshold for pruning (angle in degrees)")
                ("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
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

        string curveFilename = vm["curve"].as<std::string>();
        string inputFilename = vm["input"].as<std::string>();
        string outFilename = vm["output"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();
        double threshold = vm["thresholdPruning"].as<double>();


        GrayLevelImage2D img  = GenericReader<GrayLevelImage2D>::import( argv[ 1 ] );
        double           mass = atof( argv[ 2 ] );
        double           rmax = atof( argv[ 3 ] );
        FloatImage2D     fimg( img.domain() );
        FloatImage2D::Iterator outIt = fimg.begin();
        for ( GrayLevelImage2D::ConstIterator it = img.begin(), itE = img.end();
              it != itE; ++it )
        {
                float v = ((float)*it) / 255.0;
                *outIt++ = v;
        }
        trace.beginBlock( "Computing delta-distance." );
        Distance     delta( mass, fimg, rmax );
        const FloatImage2D& d2 = delta.myDistance2;
        trace.endBlock();

        float m = 0.0f;
        for ( typename Domain::ConstIterator it = d2.domain().begin(),
                      itE = d2.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                float v = sqrt( d2( p ) );
                m = std::max( v, m );
        }

        GradientColorMap<float> cmap_grad( 0, m );
        cmap_grad.addColor( Color( 255, 255, 255 ) );
        cmap_grad.addColor( Color( 255, 255, 0 ) );
        cmap_grad.addColor( Color( 255, 0, 0 ) );
        cmap_grad.addColor( Color( 0, 255, 0 ) );
        cmap_grad.addColor( Color( 0,   0, 255 ) );
        cmap_grad.addColor( Color( 0,   0, 0 ) );
        Board2D board;
        board << SetMode( d2.domain().className(), "Paving" );


        for ( typename Domain::ConstIterator it = d2.domain().begin(),
                      itE = d2.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                float v = sqrt( d2( p ) );
                v = std::min( (float)m, std::max( v, 0.0f ) );
                RealVector grad = delta.projection( p );
                board <<  CustomStyle( p.className(),
                                       new CustomColors( Color::Green, cmap_grad( v ) ) ) <<
                        Point(p[0] + grad[0], p[1] + grad[1]);
                board << CustomStyle( p.className(),
                                      new CustomColors( Color::Black, cmap_grad( v ) ) )
                      << p;


                // / ( 1.1 - ( (double)img( *it ) ) / 255.0 ) ;

                board.drawLine( p[ 0 ], p[ 1 ], p[ 0 ] + grad[ 0 ], p[ 1 ] + grad[ 1 ], 0 );

        }
        std::cout << endl;
        board.saveSVG("delta2.svg");
        return 0;
}
