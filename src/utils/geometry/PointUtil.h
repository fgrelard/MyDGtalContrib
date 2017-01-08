#ifndef POINT_UTIL_H
#define POINT_UTIL_H

#include <vector>
#include <set>
#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/graph/DepthFirstVisitor.h"

namespace PointUtil {
	template <typename Domain, typename Container>
	Domain computeBoundingBox(const Container& points);


	template <typename Point, typename Container, typename Vector>
	Point trackPoint(const Point& initial, const Container& container, const Vector& vector);
}


template <typename Domain, typename Container>
Domain PointUtil::computeBoundingBox(const Container & points) {
	typedef typename Container::Point Point;
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



template <typename Point, typename Container, typename Vector>
Point PointUtil::trackPoint(const Point& initial, const Container& container, const Vector& vector) {
	Point point(initial);
	int scalar = 1;
	while (container.find(point) != container.end()) {
		point = initial + vector*scalar;
		scalar++;
	}
	return point;
}

#endif
