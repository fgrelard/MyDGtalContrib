#ifndef POINT_UTIL_H
#define POINT_UTIL_H

#include <vector>
#include <set>
#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "geometry/SetProcessor.h"

namespace PointUtil {
	template <typename Domain, typename Container>
	Domain computeBoundingBox(const Container& points);

	template <typename Container, typename Point, typename RealVector>
	Point trackPoint(const Container& container, const Point& start, const RealVector& vector);

	template <typename Container, typename Point, typename RealVector>
	Container traversedLineTracking(const Container& container, const Point& start, const RealVector& vector);

	template <typename Container, typename Point, typename RealVector>
	std::pair<Point, Point> twoClosestPointsTrackingWithNormal(const Container& container,
															   const Point& ref, const RealVector& vRef,
															   const Point& other, const RealVector& vOther);
}


template <typename Domain, typename Container>
Domain PointUtil::computeBoundingBox(const Container & points) {
	typedef typename Container::value_type Point;
	int maximum = std::numeric_limits<int>::max();
	int min_x = maximum, min_y = maximum, min_z = maximum;
	int max_x = -maximum, max_y = -maximum, max_z = -maximum;

	Point low, up;
	for (const auto & point : points) {
		for (int i = 0; i < Point::dimension; i++) {
			low[i] = point[i] < low[i] ? point[i] : low[i];
			up[i] = point[i] > up[i] ? point[i] : up[i];
		}
	}
	Domain domain(low, up);
	return domain;
}

template <typename Container, typename Point, typename RealVector>
Point
PointUtil::
trackPoint(const Container& container, const Point& start, const RealVector& vector) {
	Point point(start);
	int scalar = 1;
	while (container.find(point) != container.end()) {
		point = start + vector*scalar;
		scalar++;
	}
	return point;
}

template <typename Container, typename Point, typename RealVector>
Container
PointUtil::
traversedLineTracking(const Container& volume, const Point& start, const RealVector& dir) {
	Container container(volume.domain());
	RealVector vector = dir.getNormalized();
	Point tracked = trackPoint(volume, start, vector);
	BresenhamAlgorithm<Point> bresenham(start, tracked);
	std::vector<Point> path = bresenham.linkPoints();
	container.insert(path.begin(), path.end());
	return container;
}



template <typename Container, typename Point, typename RealVector>
std::pair<Point, Point>
PointUtil::
twoClosestPointsTrackingWithNormal(const Container& container, const Point& reference, const RealVector& dirRef,
								   const Point& other, const RealVector& dirOther) {

	RealVector normalRef = dirRef.getNormalized();
	RealVector normalOther = dirOther.getNormalized();
	RealVector dirVectorReference = (other - reference).getNormalized();
	RealVector dirVectorOther = (reference - other).getNormalized();
	if (normalRef.dot(dirVectorReference) < 0)
		normalRef = -normalRef;
	if (normalOther.dot(dirVectorOther) < 0)
		normalOther = -normalOther;

	Container traversedCurrent = PointUtil::traversedLineTracking(container, reference, normalRef);
	Container traversedReference = PointUtil::traversedLineTracking(container, other, normalOther);
	double distanceCR = std::numeric_limits<double>::max();
	Point closest1, closest2;
	DGtal::ExactPredicateLpSeparableMetric<typename Container::Space,2> l2Metric;

	SetProcessor<Container> setProcessor(traversedReference);
	for (auto it = traversedCurrent.begin(), ite = traversedCurrent.end(); it != ite; ++it) {

		Point nearest = setProcessor.closestPointAt(*it);
		double currentDistance = l2Metric(nearest, *it);
		if (currentDistance < distanceCR && currentDistance > sqrt(3)) {
			distanceCR = currentDistance;
			closest1 = *it;
			closest2 = nearest;
		}
	}
	return std::make_pair(closest1, closest2);
}


#endif
