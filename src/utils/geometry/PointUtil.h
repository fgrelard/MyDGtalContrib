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
	int maximum = std::numeric_limits<int>::max();
	int min_x = maximum, min_y = maximum, min_z = maximum;
	int max_x = -maximum, max_y = -maximum, max_z = -maximum;
	for (const auto & point : points) {
		min_x = point[0] < min_x ? point[0] : min_x;
		min_y = point[1] < min_y ? point[1] : min_y;
		min_z = point[2] < min_z ? point[2] : min_z;
		max_x = point[0] > max_x ? point[0] : max_x;
		max_y = point[1] > max_y ? point[1] : max_y;
		max_z = point[2] > max_z ? point[2] : max_z;
	}
	Domain domain({min_x-1, min_y-1, min_z-1}, {max_x+1, max_y+1, max_z+1});
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
