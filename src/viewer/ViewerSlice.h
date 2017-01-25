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

public:
        ViewerSlice() = delete;
        ViewerSlice(const ViewerSlice& other) = delete;

        ViewerSlice( const std::vector<Plane>& planes,
                    const Image& reference,
                    int patchWidth);

public:
        void keyPressEvent(QKeyEvent * e);
private:
        std::vector<Plane> myPlanes;
        Image myImage;
        int myPatchWidth;
private:
        Domain myDomain3D;
        SubDomain myDomain2D;
        int myCurrentIndex = 0;

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
keyPressEvent(QKeyEvent * e) {
        typedef typename Space::RealPoint RealPoint;
        bool display = false;
        if (e->key() == Qt::Key_Right) {
                if (myCurrentIndex < myPlanes.size() - 1)
                        myCurrentIndex++;
                display = true;
        }
        if (e->key() == Qt::Key_Left) {
                if (myCurrentIndex > 0)
                        myCurrentIndex--;
                display = true;
        }
        if (display) {
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

                //Removes the previous image
                if (numberImage - 1 >= 0) {
                        (*this) << DGtal::UpdateImage3DEmbedding<Space, KSpace>(numberImage-1,
                                                                                embedder(RealPoint(0,0)),
                                                                                embedder(RealPoint(0,0)),
                                                                                embedder(RealPoint(0,0)),
                                                                                embedder(RealPoint(0,0)));
                }

                (*this) << DGtal::Viewer3D<>::updateDisplay;

                //Keeps the camera parameters (pos and zoom)
                camera->setPosition(pos);
                camera->setSceneCenter(sceneCenter);
                camera->setSceneRadius(sceneRadius);
                this->setCamera(camera);
        }

}


#endif
