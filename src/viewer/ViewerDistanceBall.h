#ifndef VIEWER_DISTANCE_BALL_H
#define VIEWER_DISTANCE_BALL_H

#include <DGtal/kernel/SpaceND.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h>
#include <DGtal/topology/KhalimskySpaceND.h>
#include <QGLViewer/qglviewer.h>
#include <shapes/Ball.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include <DGtal/io/DrawWithDisplay3DModifier.h>
#include <DGtal/io/Display3D.h>

template <typename Distance,
        typename Space = DGtal::SpaceND<3>,
        typename KSpace = DGtal::KhalimskySpaceND<3> >
class ViewerDistanceBall : public DGtal::Viewer3D<Space, KSpace> {

public:
    typedef DGtal::Viewer3D<Space, KSpace> Base;
    typedef DGtal::HyperRectDomain<Space> Domain;
    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image;
    typedef typename Space::Point Point;

    typedef DGtal::SpaceND<Space::dimension - 1> SubSpace;
    typedef DGtal::HyperRectDomain<SubSpace> SubDomain;
    typedef DGtal::ImageContainerBySTLVector<SubDomain, unsigned char> SubImage;
    typedef typename SubSpace::Point SubPoint;
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

public:
    ViewerDistanceBall() = delete;

    ViewerDistanceBall(const ViewerDistanceBall &other) = delete;

    ViewerDistanceBall(const Distance &distance);

public:

    void addVoxel(const Point &p);

    virtual void postSelection(const QPoint &point);

    virtual void mouseMoveEvent(QMouseEvent *e);

    virtual void mousePressEvent(QMouseEvent *e);

    virtual void mouseReleaseEvent(QMouseEvent *e);

    Distance getMyDistance() const;


private:
    Distance myDistance;
private:
    Domain myDomain;
    bool myRemove = true;

};

template <typename Distance, typename Space, typename KSpace>
void
ViewerDistanceBall<Distance, Space, KSpace>::
addVoxel(const Point &p) {
    Point rp = DGtal::Display3D<Space, KSpace>::embed(p);
    double width = 0.5;
    this->updateBoundingBox(rp);
    typename DGtal::Display3D<Space, KSpace>::CubeD3D v;
    v.center = rp;
    v.color = this->getFillColor();
    v.width = width;
    v.name = this->name3d();
    this->myCubesMap[v.name].insert(this->myCubesMap[v.name].begin(), v);
}


template <typename Distance, typename Space, typename KSpace>
ViewerDistanceBall<Distance, Space, KSpace>::
ViewerDistanceBall(const Distance &distance) : myDistance(distance) {
    auto lower = myDistance.domain().lowerBound();
    auto upper = myDistance.domain().upperBound();
    myDomain = Domain(Point(lower[0], lower[1], 0),
                      Point(upper[0], upper[1], 0));
}


template <typename Distance, typename Space, typename KSpace>
void
ViewerDistanceBall<Distance, Space, KSpace>::
postSelection(const QPoint &point) {
    // Compute orig and dir, used to draw a representation of the intersecting line
    auto camera = this->camera();
    auto pos = camera->position();
    auto sceneCenter = camera->sceneCenter();
    auto sceneRadius = camera->sceneRadius();
    qglviewer::Vec selectedPoint, orig, dir;
    camera->convertClickToLine(point, orig, dir);

    // Find the selectedPoint coordinates, using camera()->pointUnderPixel().
    bool found;
    L2Metric l2Metric;
    selectedPoint = camera->pointUnderPixel(point, found);
    if (found) {
        Point dgtalPoint(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
        auto iterator = std::min_element(myDomain.begin(), myDomain.end(), [&](const Point &f, const Point &s) {
            return (l2Metric(f, dgtalPoint) < l2Metric(s, dgtalPoint));
        });
        Point center = *iterator;
        SubPoint p(center[0], center[1]);
        Ball<Point> ball(center, myDistance.distance(p));
        (*this) << DGtal::CustomColors3D(DGtal::Color::Black, DGtal::Color::Black);
        for (const Point &p  : ball.pointSet()) {
            if (p[2] == 0)
                addVoxel(p);
        }

        (*this) << DGtal::Viewer3D<Space, KSpace>::updateDisplay;

        camera->setPosition(pos);
        camera->setSceneCenter(sceneCenter);
        camera->setSceneRadius(sceneRadius);
        this->setCamera(camera);
    }
}

template <typename Distance, typename Space, typename KSpace>
void
ViewerDistanceBall<Distance, Space, KSpace>::
mouseMoveEvent(QMouseEvent *e) {
    QGLViewer::mouseMoveEvent(e);
}


template <typename Distance, typename Space, typename KSpace>
void
ViewerDistanceBall<Distance, Space, KSpace>::
mousePressEvent(QMouseEvent *e) {
    QGLViewer::mousePressEvent(e);
}


template <typename Distance, typename Space, typename KSpace>
void
ViewerDistanceBall<Distance, Space, KSpace>::
mouseReleaseEvent(QMouseEvent *e) {
    QGLViewer::mouseReleaseEvent(e);

}


#endif
