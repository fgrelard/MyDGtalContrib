#ifndef ORTHOGONAL_PLANE_ESTIMATOR_H
#define ORTHOGONAL_PLANE_ESTIMATOR_H

#include "geometry/VCMAdjustableRadius.h"
#include "shapes/DigitalPlane.h"

template <typename Container, typename KernelFunction>
class OrthogonalPlaneEstimator {
private:
        static Container emptyContainer;
public:
        typedef typename Container::value_type Point;
        typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology< Adj26, Adj6 > DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;
        typedef typename Space::RealVector RealVector;
        typedef DGtal::HyperRectDomain< Space > Domain;
        typedef DGtal::ExactPredicateLpSeparableMetric<Space,2> L2Metric;
        typedef VCMAdjustableRadius<Space, L2Metric> VCM;
        typedef DigitalPlane<Space> Plane;
        typedef DGtal::EigenDecomposition<Point::dimension, double> LinearAlgebraTool;

public:
        OrthogonalPlaneEstimator() = delete;
        OrthogonalPlaneEstimator(const Container& container, const KernelFunction& chi, double R, double r, int aConnexity = 26, const L2Metric& l2Metric= L2Metric(), bool verbose = false);
        OrthogonalPlaneEstimator(const OrthogonalPlaneEstimator& other);

        ~OrthogonalPlaneEstimator();

private:
        bool isConvergenceReached(const Container& volume,
                                  const Plane& digPlane,
                                  double currentRadius,
                                  const RealVector& dirVector = RealVector::zero);

public:
        Plane planeAt(const Point& point,
                      const RealVector& dirVector = RealVector::zero,
                      const Container& points = emptyContainer);

        Plane convergentPlaneAt(const Point& point,
                                const Container& volume,
                                double maxRadius,
                                const RealVector& dirVector = RealVector::zero,
                                const Container& points = emptyContainer);

        void setRadius(double radius);

public:
        OrthogonalPlaneEstimator& operator=(const OrthogonalPlaneEstimator& other);



protected:
        VCM* myVCM;
        KernelFunction myChi;
        int myConnexity;
};

template <typename Container, typename KernelFunction>
Container OrthogonalPlaneEstimator<Container, KernelFunction>::emptyContainer = Container(Domain(Point(0,0,0), Point(0,0,0)));


template <typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
OrthogonalPlaneEstimator(const Container& container, const KernelFunction& chi, double R, double r, int aConnexity, const L2Metric& l2Metric, bool verbose) : myChi(chi), myConnexity(aConnexity) {
        myVCM = new VCM(R, r, l2Metric, verbose);
        myVCM->init(container.begin(), container.end());
}

template <typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
OrthogonalPlaneEstimator(const OrthogonalPlaneEstimator& other) {
        myVCM = new VCM(*other.myVCM);
        myChi = other.myChi;
        myConnexity = other.myConnexity;
}

template <typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
~OrthogonalPlaneEstimator() {
        if (myVCM != 0) {
                delete myVCM;
                myVCM = 0;
        }
}

template <typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>&
OrthogonalPlaneEstimator<Container, KernelFunction>::
operator=(const OrthogonalPlaneEstimator& other) {
        myVCM = new VCM( *other.myVCM );
        myChi = other.myChi;
        myConnexity = other.myConnexity;
        return *this;
}

template <typename Container, typename KernelFunction>
typename OrthogonalPlaneEstimator<Container, KernelFunction>::Plane
OrthogonalPlaneEstimator<Container, KernelFunction>::
planeAt(const Point& point, const RealVector& dirVector,const Container& points) {

        typename LinearAlgebraTool::Matrix vcm_r, evec;
        RealVector eval;
        // Compute VCM and diagonalize it.
        if (dirVector != RealVector::zero)
                vcm_r = myVCM->measureJunction( dirVector, myChi, point);
        else if (!points.empty())
                vcm_r = myVCM->measure(points, myChi, point);
        else
                vcm_r = myVCM->measure( myChi, point);
        LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
        // Display normal
        RealVector normal = evec.column(0);
        Plane plane(point, normal, myConnexity);
        return plane;
}

template <typename Container, typename KernelFunction>
bool
OrthogonalPlaneEstimator<Container, KernelFunction>::
isConvergenceReached(const Container& volume,
                     const Plane& plane,
                     double currentRadius,
                     const RealVector& dirVector) {

        L2Metric l2Metric;
        bool alright = true;
        typename Plane::DigitalSet intersection = plane.intersectionWithSetOneCC(volume);
        if (dirVector == RealVector::zero) {
                for (const Point& p : intersection) {
                        if (l2Metric(p, plane.getCenter()) >= currentRadius) {
                                alright = false;
                        }
                }
        } else {
                Point current = plane.getCenter();
                int scalar = 1;
                auto itSetVolume = std::find(volume.begin(), volume.end(), current);
                while (itSetVolume != volume.end()) {
                        current = plane.getCenter() + dirVector * scalar;
                        scalar++;
                        itSetVolume = std::find(volume.begin(), volume.end(), current);
                        double distance = l2Metric(current, plane.getCenter());
                        alright = (distance < currentRadius);
                }
                return alright;
        }
}

template <typename Container, typename KernelFunction>
typename OrthogonalPlaneEstimator<Container, KernelFunction>::Plane
OrthogonalPlaneEstimator<Container, KernelFunction>::
convergentPlaneAt(const Point& point,
                  const Container& volume,
                  double maxRadius,
                  const RealVector& dirVector,
                  const Container& points) {

        bool isConvergent = false;
        double currentRadius = myVCM->r();
        Plane convergentPlane;
        do {
                currentRadius++;
                setRadius(currentRadius);
                convergentPlane = planeAt(point, dirVector, points);
                isConvergent = isConvergenceReached(volume,
                                                    convergentPlane,
                                                    currentRadius,
                                                    dirVector);

        } while (!isConvergent && currentRadius < maxRadius);

        return convergentPlane;
}

template <typename Container, typename KernelFunction>
void
OrthogonalPlaneEstimator<Container, KernelFunction>::
setRadius(double radius) {
        myVCM->setMySmallR(radius);
        myChi = KernelFunction(1.0, radius);
}

#endif
