#ifndef ORTHOGONAL_PLANE_ESTIMATOR_H
#define ORTHOGONAL_PLANE_ESTIMATOR_H

#include "geometry/VCMAdjustableRadius.h"
#include "shapes/DigitalPlane.h"

template <typename Container, typename KernelFunction>
class OrthogonalPlaneEstimator {
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
        OrthogonalPlaneEstimator(const Container& container, const KernelFunction& chi, double R, double r, const L2Metric& l2Metric= L2Metric(), bool verbose = false);
        OrthogonalPlaneEstimator(const OrthogonalPlaneEstimator& other);

        ~OrthogonalPlaneEstimator();

public:
        Plane planeAt(const Point& point,
                      const RealVector& dirVector = RealVector::zero,
                      const Container& points = Container(Domain(Point(0,0,0), Point(0,0,0)));

public:
        OrthogonalPlaneEstimator& operator=(const OrthogonalPlaneEstimator& other);



protected:
        VCM* myVCM;
        KernelFunction myChi;
};


template <typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
OrthogonalPlaneEstimator(const Container& container, const KernelFunction& chi, double R, double r, const L2Metric& l2Metric, bool verbose) : myChi(chi) {
        myVCM = new VCM(R, r, l2Metric, verbose);
        myVCM->init(container.begin(), container.end());
}

template <typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
OrthogonalPlaneEstimator(const OrthogonalPlaneEstimator& other) {
        myVCM = new VCM(*other.myVCM);
        myChi = other.myChi;
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
                vcm_r = myVCM.measureJunction( dirVector, myChi, point);
        else if (!points.empty())
                vcm_r = myVCM.measure(points, myChi, point);
        else
                vcm_r = myVCM.measure( myChi, point);
        LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
        // Display normal
        RealVector normal = evec.column(0);
        Plane plane(point, normal);
        return plane;
}


#endif
