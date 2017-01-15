#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/MultiPathThinner.h"

using namespace DGtal;
using namespace std;


void testMultiThin(Viewer3D<>& viewer) {
        using namespace Z3i;
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;

        std::string home = getenv("HOME");
        std::string path = home + "/test_img/cylinder.vol";
        Image image = GenericReader<Image>::import(path);
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        Point ref(0,0,25);
        RealVector normal(0.1,0.2,-0.9);

        Point toLink1(3,3,3);
        RealVector normal2(-0.1,0,0.9);

        Point toLink2(-2,-2,-2);
        RealVector normal3(-0.1,0,0.9);

        Point toLink3(0, -4, 10);
        RealVector normal4(0, 1,0);

        Point toLink4(0, 4, 10);
        RealVector normal5(0, -1,0);

        pair<Point, RealVector> pairReference = std::make_pair(ref, normal);
        pair<Point, RealVector> pairToLink = std::make_pair(toLink1, normal2);
        pair<Point, RealVector> pairToLink2 = std::make_pair(toLink2, normal3);
        pair<Point, RealVector> pairToLink3 = std::make_pair(toLink3, normal4);
        pair<Point, RealVector> pairToLink4 = std::make_pair(toLink4, normal5);
        std::vector<std::pair<Point, RealVector> > vectors;
        vectors.push_back(pairToLink);
        vectors.push_back(pairToLink2);
        vectors.push_back(pairToLink3);
        vectors.push_back(pairToLink4);
        MultiPathThinner<Z3i::DigitalSet> multi(setVolume, vectors, pairReference);
        viewer << CustomColors3D(Color::Red, Color::Red) << ref << toLink1 << toLink2<<toLink3<<toLink4;
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << multi.linkPointsThin();
        viewer << Viewer3D<>::updateDisplay;

}


int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testMultiThin(viewer);
        app.exec();
        return 0;
}
