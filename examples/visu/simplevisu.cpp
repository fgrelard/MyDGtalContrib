#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "shapes/DigitalPlane.h"
#include "geometry/DigitalPlaneProcessor.h"


using namespace std;
using namespace DGtal;



int main( int  argc, char**  argv )
{

    QApplication application(argc,argv);
    Viewer3D<> viewer;
    viewer.show();

    viewer << CustomColors3D(Color::Blue, Color::Blue) << Z3i::Point(0,0,0) << Z3i::Point(1,1,1);
    Z3i::RealVector n(1,1,1);
     Z3i::Point current(1,1,1);
     DigitalPlane<Z3i::Space> digPlane(current, n.getNormalized());
    DigitalPlaneProcessor<Z3i::Space> digProc(digPlane);
    std::vector<Z3i::RealVector> plane = digProc.planeToQuadrangle();
    double f = 20;

    Z3i::RealPoint current2(1,1,1);
    viewer << CustomColors3D(Color::Red, Color::Red);
    viewer.addQuad(current2 + plane[0] * f,
                   current2 + plane[1] * f,
                   current2 + plane[2] * f,
                   current2 + plane[3] * f);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j <  3; j++) {
            for (int k = 0; k < 3; k++) {
                Color c(0,0,255,30);
                Z3i::Point p(i,j,k);
                if (digPlane.isPointAbove(p))
                    viewer << CustomColors3D(Color::Yellow, Color::Yellow) << p;
                else
                    viewer << CustomColors3D(c,c) << p;
            }
        }
    }

    viewer << Viewer3D<>::updateDisplay;
    application.exec();
    return 0;


}
