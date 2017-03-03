#include <iostream>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "shapes/Ball.h"
#include "vcm/delta/SegmentationAnnulus.h"


using namespace std;
using namespace DGtal;


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

        Ball<Point> ball(Point(40,60), 10);
        DigitalSet set = ball.pointSet();
        SegmentationAnnulus<FloatImage2D> segm(delta, set);
        DigitalSet annulus = segm.extractAnnulus();
        GrayLevelImage2D outImage(img.domain());
        imageFromRangeAndValue(annulus.begin(), annulus.end(), outImage, 255);
        PGMWriter<GrayLevelImage2D>::exportPGM("annulus.pgm", outImage);
        return 0;
}
