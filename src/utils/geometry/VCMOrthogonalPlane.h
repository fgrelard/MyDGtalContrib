#ifndef VCM_ORTHOGONAL_PLANE_H
#define VCM_ORTHOGONAL_PLANE_H

#include "DGtal/geometry/surfaces/ParallelStrip.h"
#include "VoronoiCovarianceMeasure.h"

template <typename TSpace>
class VCMOrthogonalPlane {
public:
        typedef TSpace Space;
        typedef typename Space::Point Point;
public:
        VCMOrthogonalPlane() {}
        template <typename VCM, typename KernelFunction>
    DGtal::Z3i::RealPoint computeNormalFromVCM(const DGtal::Z3i::Point& currentPoint, const VCM& vcm, const KernelFunction& chi,
                                               int coordinate, const DGtal::Z3i::RealVector& dirVector = DGtal::Z3i::RealVector(),
                                               const std::vector<DGtal::Z3i::Point>& v = std::vector<DGtal::Z3i::Point>());

private:
        DGtal::ParallelStrip<TSpace> myPlane;
        Point myPlaneCenter;
};

#endif
