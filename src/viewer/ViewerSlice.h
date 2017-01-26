#ifndef VIEWER_SLICE_H
#define VIEWER_SLICE_H

#include "DGtal/kernel/SpaceND.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/viewers/DrawWithViewer3DModifier.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/base/BasicFunctors.h"
#include "shapes/DigitalPlane.h"
#include "geometry/DigitalPlaneProcessor.h"


template <typename Space = DGtal::SpaceND<3>,
          typename KSpace = DGtal::KhalimskySpaceND<3> >
class ViewerSlice : public DGtal::Viewer3D<Space, KSpace> {
public:
        typedef DGtal::HyperRectDomain<Space> Domain;
        typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image;
        typedef DigitalPlane<Space> Plane;
        typedef typename Space::Point Point;

        typedef DGtal::SpaceND<Space::dimension-1> SubSpace;
        typedef DGtal::HyperRectDomain<SubSpace> SubDomain;
        typedef DGtal::ImageContainerBySTLVector<SubDomain, unsigned char> SubImage;
        typedef typename SubSpace::Point SubPoint;
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

public:
        ViewerSlice() = delete;
        ViewerSlice(const ViewerSlice& other) = delete;

        ViewerSlice( const std::vector<Plane>& planes,
                    const Image& reference,
                    int patchWidth);

public:
        virtual void keyPressEvent(QKeyEvent * e);
        virtual void postSelection(const QPoint& point);
        virtual void mouseMoveEvent ( QMouseEvent *e );
        virtual void mousePressEvent ( QMouseEvent *e );
        virtual void mouseReleaseEvent ( QMouseEvent *e );

private:
        void displayCurrentSlice();
private:
        std::vector<Plane> myPlanes;
        Image myImage;
        int myPatchWidth;
private:
        Domain myDomain3D;
        SubDomain myDomain2D;
        int myCurrentIndex = 0;
        bool myRemove = true;

};


template <typename Space, typename KSpace>
ViewerSlice<Space, KSpace>::
ViewerSlice(const std::vector<Plane>& planes,
            const Image& reference,
            int patchWidth) : myPlanes(planes), myPatchWidth(patchWidth), myImage(reference) {
        Point upper = reference.domain().upperBound();
        Point lower = reference.domain().lowerBound();
        myDomain3D = Domain(lower - Point::diagonal(myPatchWidth),
                            upper + Point::diagonal(myPatchWidth));

        myDomain2D = SubDomain(SubPoint::zero,
                               SubPoint::diagonal(myPatchWidth));
}

template <typename Space, typename KSpace>
void
ViewerSlice<Space, KSpace>::
displayCurrentSlice() {
        typedef typename Space::RealPoint RealPoint;

        auto camera = this->camera();
        auto pos = camera->position();
        auto sceneCenter = camera->sceneCenter();
        auto sceneRadius = camera->sceneRadius();
        Plane plane = myPlanes[myCurrentIndex];
        DGtal::functors::Identity idV;
        DGtal::functors::Point2DEmbedderIn3D<Domain>  embedder(myDomain3D,
                                                               plane.getCenter(),
                                                               plane.getPlaneEquation().normal(),
                                                               myPatchWidth,
                                                               myDomain3D.lowerBound());
        SubImage image = DigitalPlaneProcessor<Space>(plane).sliceFromPlane(myImage, myPatchWidth);
        (*this) << image;
        int numberImage =  this->getCurrentGLImageNumber() - 1;

        (*this) << DGtal::UpdateImage3DEmbedding<Space, KSpace>(numberImage,
                                                                embedder(RealPoint(0,0)),
                                                                embedder(RealPoint(myPatchWidth,0)),
                                                                embedder(myDomain2D.upperBound()),
                                                                embedder(RealPoint(0, myPatchWidth)));

        //Removes the previous images
        if (myRemove) {
                for (int i = 0; i < numberImage; i++) {
                        (*this) << DGtal::UpdateImage3DEmbedding<Space, KSpace>(i,
                                                                                embedder(RealPoint(0,0)),
                                                                                embedder(RealPoint(0,0)),
                                                                                embedder(RealPoint(0,0)),
                                                                                embedder(RealPoint(0,0)));
                }
        }

        (*this) << DGtal::Viewer3D<>::updateDisplay;

        //Keeps the camera parameters (pos and zoom)
        camera->setPosition(pos);
        camera->setSceneCenter(sceneCenter);
        camera->setSceneRadius(sceneRadius);
        this->setCamera(camera);
}

template <typename Space, typename KSpace>
void
ViewerSlice<Space, KSpace>::
keyPressEvent(QKeyEvent * e) {

        bool display = false;
        if (e->key() == Qt::Key_Right) {
                if (myCurrentIndex < myPlanes.size() - 1)
                        myCurrentIndex++;
                display = true;
        }
        else if (e->key() == Qt::Key_Left) {
                if (myCurrentIndex > 0)
                        myCurrentIndex--;
                display = true;
        }
        else if (e->key() == Qt::Key_R) {
                myRemove = !myRemove;
        }
        else {
                DGtal::Viewer3D<>::keyPressEvent(e);
        }
        if (display) {
                displayCurrentSlice();
        }

}

template <typename Space, typename KSpace>
void
ViewerSlice<Space,KSpace>::
postSelection(const QPoint& point)
{
  // Compute orig and dir, used to draw a representation of the intersecting line
        auto camera = this->camera();
        qglviewer::Vec selectedPoint, orig, dir;
        camera->convertClickToLine(point, orig, dir);

        // Find the selectedPoint coordinates, using camera()->pointUnderPixel().
        bool found;
        L2Metric l2Metric;
        selectedPoint = camera->pointUnderPixel(point, found);
        if (found) {
                Point dgtalPoint(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
                auto iterator = std::min_element(myPlanes.begin(), myPlanes.end(), [&](const Plane& f, const Plane& s) {
                                return (l2Metric(f.getCenter(), dgtalPoint) < l2Metric(s.getCenter(), dgtalPoint));
                        });
                myCurrentIndex = iterator - myPlanes.begin();
                displayCurrentSlice();
        }
}

template <typename Space, typename KSpace>
void
ViewerSlice<Space, KSpace>::
mouseMoveEvent ( QMouseEvent *e ) {
        QGLViewer::mouseMoveEvent(e);
}


template <typename Space, typename KSpace>
void
ViewerSlice<Space, KSpace>::
mousePressEvent ( QMouseEvent *e ) {
        QGLViewer::mousePressEvent(e);
}

template <typename Space, typename KSpace>
void
ViewerSlice<Space, KSpace>::
mouseReleaseEvent ( QMouseEvent *e ) {
        QGLViewer::mouseReleaseEvent(e);

}


#endif
