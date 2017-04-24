#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/viewers/Viewer3D.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int argc, char** argv )
{
        using namespace DGtal;
        using namespace DGtal::Z2i;

        typedef ImageContainerBySTLVector<Domain,unsigned char> GrayLevelImage2D;
        typedef ImageContainerBySTLVector<Domain,float>         FloatImage2D;
    typedef DistanceToMeasureEdge <FloatImage2D> Distance;

        GrayLevelImage2D img  = GenericReader<GrayLevelImage2D>::import( "/home/florent/Documents/DGtal/VCM/img/annulus.pgm" );
        double           mass = 3;
        double           rmax = 10;
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
        trace.endBlock();

        float m = 0.0f;
        for ( typename Domain::ConstIterator it = delta.domain().begin(),
                      itE = delta.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                float v = sqrt( delta.distance2( p ) );
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
        board << SetMode( delta.domain().className(), "Paving" );


        for ( typename Domain::ConstIterator it = delta.domain().begin(),
                      itE = delta.domain().end(); it != itE; ++it )
        {
                Point p = *it;
                float v = sqrt( delta.distance2( p ) );
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
