#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "geometry/path/LinkPointAlgorithm.h"
#include "geometry/path/BezierCasteljauLinkAlgorithm.h"
#include "geometry/path/BezierLinkAlgorithm.h"
#include "geometry/path/BresenhamAlgorithm.h"
#include "geometry/path/DijkstraAlgorithm.h"
#include "geometry/path/AStarAlgorithm.h"
#include "shapes/Border.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/SetProcessor.h"

using namespace DGtal;
using namespace std;

void testMajorAxis() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        SetProcessor<Z3i::DigitalSet> setProcessor(setVolume);
        pair<Z3i::Point, Z3i::Point> pair = setProcessor.majorAxis();
        trace.info() << pair.first << " " << pair.second << endl;
}


void testLengthMajorAxis() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;


        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        SetProcessor<Z3i::DigitalSet> setProcessor(setVolume);
        double radius = setProcessor.lengthMajorAxis();
        trace.info() << radius << endl;
        const Color CURVE3D_COLOR( 100, 100, 140, 128 );
}


void testClosestPointAt() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;


        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        SetProcessor<Z3i::DigitalSet> setProcessor(setVolume);
        Z3i::RealPoint realPoint(1.2, 2.8, 32.5);
        Z3i::Point closest = setProcessor.closestPointAt(realPoint);
        trace.info() << closest << endl;
        const Color CURVE3D_COLOR( 100, 100, 140, 128 );
}


void testSameContainer() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;


        Image image = GenericReader<Image>::import("/home/florent/test_img/thskeleton_boudin.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        SetProcessor<Z3i::DigitalSet> setProcessor(setVolume);
        bool isSame1 = setProcessor.sameContainer(setVolume);
        setVolume.erase(setVolume.begin());
        bool isSame2 = setProcessor.sameContainer(setVolume);
        trace.info() << "Test " << ((isSame1 && !isSame2) ? " passed" : " failed") << endl;
}

int main(int argc, char** argv) {
        testMajorAxis();
        testLengthMajorAxis();
        testClosestPointAt();
        testSameContainer();
        return 0;
}
