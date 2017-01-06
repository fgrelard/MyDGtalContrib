#ifndef POINT_UTIL_H
#define POINT_UTIL_H

#include <vector>
#include <set>
#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DepthFirstVisitor.h"

namespace PointUtil {
	template <typename Point>
	bool areNeighbors(const Point& point, const Point& other);


	template <typename Domain, typename Container>
	Domain computeBoundingBox(const Container& points);


	/*
	 * We have to visit the direct neighbours in order to have a container with voxels
	 * ordered sequentially by their connexity
	 * Otherwise we have a point container with points which are not neighbours
	 * and this impairs maximal segment recognition
	 */
	template <typename Container, typename Image, typename Point>
	Container containerFromDepthTraversal(const Image& image, const Point& point, int thresholdMin,
		int thresholdMax);

	template <typename Container, typename OtherContainer, typename Point>
	Container containerFromDepthTraversal(const OtherContainer& image, const Point& point);


	template <typename Point, typename Container, typename Vector>
	Point trackPoint(const Point& initial, const Container& container, const Vector& vector);
}

template <typename Point>
bool PointUtil::areNeighbors(const Point& point, const Point& other) {
	typename Point::Scalar otherx = other[0];
	typename Point::Scalar othery = other[1];
	typename Point::Scalar otherz = other[2];

	typename Point::Scalar pointx = point[0];
	typename Point::Scalar pointy = point[1];
	typename Point::Scalar pointz = point[2];

	bool sameX = pointx == otherx || pointx == otherx + 1 || pointx == otherx-1;
	bool sameY = pointy == othery || pointy == othery + 1 || pointy == othery-1;
	bool sameZ = pointz == otherz || pointz == otherz + 1 || pointz == otherz-1;
	return sameX && sameY && sameZ;
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

template <typename Container, typename Image, typename Point>
Container PointUtil::containerFromDepthTraversal(const Image& image, const Point& point, int thresholdMin,
										  int thresholdMax) {
	using namespace DGtal;


	typedef MetricAdjacency<Z3i::Space, 3> Graph;
	typedef DepthFirstVisitor<Graph, std::set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;

	Graph graph;
	Visitor visitor( graph, point );
	MyNode node;
	Container container;

	while ( !visitor.finished() )
	{
		node = visitor.current();
		if ( image.domain().isInside(node.first) &&
			 image(node.first) >= thresholdMin &&
			 image(node.first) <= thresholdMax ) { //is inside domain
		    container.push_back(node.first);
			visitor.expand();
		}
		else
			visitor.ignore();
	}
	return container;
}

template <typename Container, typename OtherContainer, typename Point>
Container PointUtil::containerFromDepthTraversal(const OtherContainer& image, const Point& point) {
	using namespace DGtal;


	typedef DGtal::Z3i::Object26_6 Graph;
	typedef DepthFirstVisitor<Graph, std::set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;

	Graph graph(DGtal::Z3i::dt26_6, image);
	Visitor visitor( graph, point );
	MyNode node;
	Container container;

	while ( !visitor.finished() )
	{
		node = visitor.current();
		container.push_back(node.first);
		visitor.expand();
	}
	return container;
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
