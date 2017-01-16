#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <ctime>
#include <vector>
#include <cstdlib>
#include "vcm/skeleton/JunctionProcessingSkeleton.h"
#include "geometry/AbovePlanePredicate.h"

using namespace DGtal;
using namespace std;

void testConstructJunctionProcessing() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        typedef AbovePlanePredicate<Z3i::Space> Predicate;
        typedef DigitalPlane<Z3i::Space> Plane;

        Image image = GenericReader<Image>::import("/home/florent/test_img/airway.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Image curve = GenericReader<Image>::import("/home/florent/test_img/skeletonVCMAirway.vol");
        Z3i::DigitalSet setCurve(curve.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setCurve, curve, 0, 255);
        std::vector<Plane> f = {Plane(*setCurve.begin(), Z3i::RealVector(0,0,1)),
                                Plane(*(++setCurve.begin()), Z3i::RealVector(0,0,-1))};

        JunctionProcessingSkeleton<Z3i::DigitalSet, Predicate> junctionProc(setCurve, setCurve, setVolume, f);
        Z3i::DigitalSet junctionProcessed = junctionProc.postProcess();
}


int main() {
        testConstructJunctionProcessing();
        return 0;
}
