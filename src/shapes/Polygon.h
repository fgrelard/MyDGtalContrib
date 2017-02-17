#ifndef POLYGON_H
#define POLYGON_H

#include <initializer_list>
#include <vector>
#include <cstddef>

template<typename Point>
class Polygon {

public:
    Polygon() {}

    template <typename Iterator>
    Polygon(Iterator it, Iterator ite) {
        myPoints = std::vector<Point>(it, ite);
    }

    Polygon(std::initializer_list<Point> l) {
        myPoints = l;
    }

    std::vector<Point> getPolygon() const { return myPoints; }

    bool isInside(const Point &p);

private:
    std::vector<Point> myPoints;
};

template <typename Point>
bool Polygon<Point>::
isInside(const Point &p) {
    bool inside = false;
    for (size_t c= 0, d = myPoints.size()-1; c < myPoints.size(); d = c++)
    {
        Point current = myPoints[c];
        Point other = myPoints[d];
        Point vector = other - current;
        if( ((current[1] > p[1]) != (other[1] > p[1])) &&
            (p[0] < (vector[0]) * (p[1] - current[1]) /
                   (vector[1]) + current[0]) )
            inside = !inside;
    }
    return inside;
}


#endif
