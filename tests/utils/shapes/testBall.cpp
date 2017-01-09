#include "shapes/Ball.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/viewers/Viewer3D.h"

using namespace DGtal;
using namespace Z3i;

void testBallContains() {
        Z3i::Point center(0,5,0);
        Ball<Z3i::Point> ball(center, 10);
        bool contains1 = ball.contains(Point(0,5,0));
        bool contains2 = ball.contains(Point (0,15,0));
        bool contains3 = ball.contains(Point(0, 16,0));
        trace.info() << "first test=" << (contains1==1) << endl;
        trace.info() << "second test=" << (contains2==1) << endl;
        trace.info() << "third test=" << (contains3==0) << endl;
}

void testBallPoints(Viewer3D<>& viewer) {
        Z3i::Point center(0,5,0);
        Ball<Z3i::Point> ball(center, 10);
        Z3i::DigitalSet points = ball.pointSet();
        viewer << CustomColors3D(Color::Red, Color::Red) << points;
}

void testBallHalfPoints(Viewer3D<>& viewer) {
        Z3i::Point center(0, 26,0);
        Ball<Z3i::Point> ball(center, 10);
        Z3i::DigitalSet points = ball.pointsInHalfBall(Z3i::RealVector(0,0,1));
        viewer << CustomColors3D(Color::Red, Color::Red) << points;
}



void testBallSurface(Viewer3D<>& viewer) {
        Z3i::Point center(0, 47 ,0);
        Ball<Z3i::Point> ball(center, 10);
        Z3i::DigitalSet points = ball.pointsSurfaceBall();
        viewer << CustomColors3D(Color::Red, Color::Red) << points;
}

void testBallIntersection(Viewer3D<>& viewer) {
        Z3i::Point centerI(5,70,5);
        Ball<Z3i::Point> ballI(centerI, 10);
        Z3i::DigitalSet setI = ballI.pointSet();
        Z3i::Point center(0, 68 ,0);
        Ball<Z3i::Point> ball(center, 10);
        Z3i::DigitalSet points = ball.intersection(setI);
        viewer << CustomColors3D(Color::Red, Color::Red) << points;
}

int main(int argc, char** argv) {
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();
        testBallContains();
        testBallPoints(viewer);
        testBallHalfPoints(viewer);
        testBallSurface(viewer);
        testBallIntersection(viewer);
        viewer << Viewer3D<>::updateDisplay;
        app.exec();
        return 0;
}
