#ifndef MODELLER_H
#define MODELLER_H

#include <vector>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "geometry/PointUtil.h"
#include "geometry/Distance.h"
#include "shapes/Ball.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"

template <typename Container>
class Modeller {

        BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

public:
        typedef typename Container::Point Point;
        typedef typename Container::Space Space;
        typedef DGtal::HyperRectDomain<Space> Domain;

public:
        Container drawCircle(float radius, float cx, float cy, float z, float increment);

        Eigen::Matrix<double, Eigen::Dynamic, 4> drawDisk(double radius, double cx, double cy, double cz, long int &row, float increment);

        Container drawEllipse(float a, float b, float cx, float cy, float cz, float increment);

        Container drawDisk(float radius, float cx, float cy, float cz, float increment);

        Container drawEllipsoid(float a, float b, float c, float cx, float cy, float cz, float increment);

        Container drawCone(int length, float radius, float increment);

        Eigen::Matrix<double, Eigen::Dynamic, 4> drawCylinder(int length, int radius,  float rotationAngle, float increment);

        Container drawCylinder(int length, float radius, float increment);

        Container drawDeformedCylinder(int length, int radius, float increment);

        Container createHelixCurve(int range, int radiusWinding, int radiusSpiral, int pitch, float increment);

        Container createHelixCurve( int range, int radius, int pitch, float increment);

        Container createStraightLine(int range, float increment = 0.01);

        void createSyntheticAirwayTree(Container & c, int branchNumber, int lengthTrachea, int z, float rotationAngle, Point firstPoint, float increment);

        Container createLogarithmicCurve(int range, float increment);

        Container createVolumeFromCurve(const Container & curve,  int ballRadius);

        Container createHalfVolumeFromCurve(const Container & curve, int ballRadius);

        Container createRotatedVolumeFromCurve(const Container & curve, int ballRadius, double angle, const Eigen::Vector3d& vector = Eigen::Vector3d::UnitX());

        Container addNoise(const Container& set, float noise);

        Container create2DCurve();

        template <typename Board>
        void create2DNaiveTangentsForVisu(const Container & points, Board& board);
};

template <typename Container>
Container Modeller<Container>::drawCircle(float radius, float cx, float cy, float z, float increment) {
        float x, y;
        double angle = 0.0;
        std::set<Point> set;
        while (angle <= 2 * M_PI){
                x = radius * cos( angle );
                y = radius * sin( angle );
                set.insert(Point((int)x+cx, (int)y+cy, (int)z));
                angle += increment;
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::drawEllipse(float a, float b, float cx, float cy, float cz, float increment) {
        std::set<Point> set;
        for (float x = cx -x - a, xend = cx + a + 1; x < xend; x+=increment) {
                for (float y = cy - b, yend = cy + b +  1; y < yend; y+=increment) {
                        if (Distance::euclideanDistance(x/a, y/b, cx/a, cy/b) <= 1) {
                                set.insert(Point(x, y, cz));
                        }
                }
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::drawEllipsoid(float a, float b, float c, float cx, float cy, float cz, float increment) {
        std::set<Point> set;
        Point center(cx/a, cy/b, cz/c);
        for (float x = cx - a, xend = cx + a + 1; x < xend; x+=increment) {
                for (float y = cy - b, yend = cy + b +  1; y < yend; y+=increment) {
                        for (float z = cz - c, zend = cz + c +  1; z < zend; z+=increment) {
                                Point current(x/a, y/b, z/c);
                                if (Distance::euclideanDistance(center, current) <= 1) {
                                        set.insert(Point(x,y,z));
                                }
                        }
                }
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Eigen::Matrix<double, Eigen::Dynamic, 4> Modeller<Container>::drawDisk(double radius, double cx, double cy, double cz, long int &row, float increment) {
        Eigen::Matrix<double, Eigen::Dynamic, 4> m;
        for (float x = cx - radius, xend = cx + radius + 1; x < xend; x+=increment) {
                for (float y = cy - radius, yend = cy + radius +  1; y < yend; y+=increment) {
                        if (Distance::euclideanDistance(x, y, cx, cy) <= radius) {
                                m(row, 0) = x;
                                m(row, 1) = y;
                                m(row, 2) = cz;
                                m(row, 3) = 1;
                                row++;
                                m.conservativeResize(row+1, 4);
                        }
                }
        }
        return m;
}

template <typename Container>
Container Modeller<Container>::drawDisk(float radius, float cx, float cy, float cz, float increment) {
        std::set<Point> set;
        for (float x = cx - radius, xend = cx + radius + 1; x < xend; x+=increment) {
                for (float y = cy - radius, yend = cy + radius +  1; y < yend; y+=increment) {
                        if (Distance::euclideanDistance(x, y, cx, cy) <= radius) {
                                set.insert(Point(x, y, cz));
                        }
                }
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::drawCone(int length, float radius, float increment) {
        std::set<Point> set;
        for (int i = 0; i < length; i++) {
                Container aSet = drawDisk(radius, 0, 0, i, increment);
                set.insert(aSet.begin(), aSet.end());
                radius -= increment;
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::drawCylinder(int length, float radius, float increment) {
        std::set<Point> set;
        for (int i = 0; i < length; i++) {
                Container aSet = drawDisk(radius, 0, 0, i, increment);
                set.insert(aSet.begin(), aSet.end());
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}



template <typename Container>
Eigen::Matrix<double, Eigen::Dynamic, 4> Modeller<Container>::drawCylinder(int length, int radius, float rotationAngle, float increment) {
        int progress = 0;
        long int row = 0;
        Eigen::Matrix<double, Eigen::Dynamic, 4> mm;
        while (progress < length) {
                if (progress % (length/10) != 0) {
                        Eigen::Matrix<double, Eigen::Dynamic, 4> m = drawDisk(radius, 0, 0, progress, row, increment);
                        mm << m;
                } else {
                        if (radius > 0.6 * length / 6) {
                                radius--;
                        }
                        Eigen::Matrix<double, Eigen::Dynamic, 4> m = drawDisk(radius, 0, 0, progress, row, increment);
                        mm << m;
                }
                progress++;
        }
        //Deletes last row
        mm.conservativeResize(row, 4);
        //Rotation with rotationAngle
        Eigen::Affine3d rot(Eigen::AngleAxisd(rotationAngle, Eigen::Vector3d::UnitX()));
        if (rotationAngle != 0.0) {
                mm = mm * rot.matrix();
        }
        return mm;
}

template <typename Container>
Container Modeller<Container>::drawDeformedCylinder(int length, int radius, float increment) {
        int a = radius;
        int b = radius;
        std::set<Point> set;
        for (int i = 0; i < length; i++) {
                Container aSet = drawEllipse(a, b, 0, 0, i, increment);
                set.insert(aSet.begin(), aSet.end());
                if (i % 3 == 0) {
                        a -= rand() % 3 - 1;
                        if (a <= 1) a = 2;
                        b += rand() % 3 - 1;
                        if (b<=1) b = 2;
                }
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}


/*
 * Radius  Winding : radius of the loop
 * RadiusSpiral: radius of the volume
 */
template <typename Container>
Container Modeller<Container>::createHelixCurve(int range, int radiusWinding, int radiusSpiral, int pitch, float increment) {
        std::set<Point> set;
        for (float i=0.0f; i < range; i+=increment) {
                float centerx = radiusWinding * cos(i/pitch);
                float centery = radiusWinding * sin(i/pitch);
                Container aSet = drawCircle(radiusSpiral, centerx, centery, i, increment);
                set.insert(aSet.begin(), aSet.end());
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}


/*
 * This function does not allow to set the spiral radius
 */
template <typename Container>
Container Modeller<Container>::createHelixCurve( int range, int radius, int pitch, float increment) {
        std::set<Point> set;
        for (float i=0.0; i < range; i+=increment) {
                int x = radius * cos(i);
                int y = radius * sin(i);
                int z = pitch * i;
                Point p(x,y,z);
                set.insert(p);
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::createStraightLine(int range, float increment) {
        std::set<Point> set;
        for (float i = 0; i < range; i+=increment) {
                Point p (0, 0, i);
                set.insert(p);
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}


template <typename Container>
void Modeller<Container>::createSyntheticAirwayTree(Container& c, int branchNumber, int lengthTrachea, int z, float rotationAngle, Point firstPoint, float increment) {
        if (branchNumber == 0) return c;
        int radius = lengthTrachea / 6;
        Eigen::Matrix<double, Eigen::Dynamic, 4> matrix(1, 4);
        std::set<Point> set;
        //Creates a rotated cylinder
        matrix = drawCylinder(lengthTrachea, radius, rotationAngle, increment);
        //Filling the vector with points from matrix
        for (int i = 0; i < matrix.innerSize(); i++) {
                if (matrix(i, 0) != 0 || matrix(i, 1) != 0 || matrix(i, 2) != 0) {
                        int posx = matrix(i, 0)+firstPoint[0];
                        //Translation after rotation
                        int posy = matrix(i, 1)+firstPoint[1];
                        int posz = matrix(i, 2)+firstPoint[2];
                        c.insert(Point(posx, posy, posz));
                }
        }
        z += lengthTrachea;
        //Determining initial starting point for new branches
        firstPoint += Point(0, lengthTrachea * sin(rotationAngle), lengthTrachea * cos(rotationAngle));
        matrix.resize(0,0);
        createSyntheticAirwayTree(branchNumber - 1, lengthTrachea * 0.7, z, rotationAngle + 0.2 * M_PI, firstPoint, increment);
        createSyntheticAirwayTree(branchNumber - 1, lengthTrachea * 0.7, z, rotationAngle - 0.2 * M_PI, firstPoint, increment);
}

template <typename Container>
Container Modeller<Container>::createLogarithmicCurve(int range, float increment) {
        std::set<Point> set;
        for (float i = 1; i < range; i+=increment) {
                float x = i;
                float y = i;
                float z = 20*log(i);
                set.insert(Point((int)x, (int)y, (int)z));
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}


template <typename Container>
Container Modeller<Container>::createVolumeFromCurve(const Container & curve, int ballRadius) {
        std::set<Ball<Point> > set;
        for (const Point& point : curve) {
                set.insert(Ball<Point>(point, ballRadius));
        }
        std::set<Point> setP;
        for (const Ball<Point>& current : set) {
                for (const Point& point : current.pointsInBall()) {
                        setP.insert(point);
                }
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::createHalfVolumeFromCurve(const Container & curve, int ballRadius) {
        std::set<Ball<Point> > set;
        for (const Point& point : curve) {
                set.insert(Ball<Point>(point, ballRadius));
        }
        std::set<Point> setP;
        for (const Ball<Point>& current : set) {
                for (const Point& point : current.pointsInHalfBall()) {
                        setP.insert(point);
                }
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}


template <typename Container>
Container Modeller<Container>::createRotatedVolumeFromCurve(const Container& curve, int ballRadius, double angle, const Eigen::Vector3d& vectorRot) {
        typedef Ball<Point> Ball;
        typedef typename Point::Scalar Scalar;
        Eigen::Matrix<double, Eigen::Dynamic, 4> matrix(1, 4);
        long int row = 0;

        //drawDisk(matrix, 20.0, 0, 0, 0, row, 0.5);
        for (const Point& point : curve) {
                matrix.conservativeResize(row+1, 4);
                matrix(row, 0) = point[0];
                matrix(row, 1) = point[1];
                matrix(row, 2) = point[2];
                matrix(row, 3) = 1;
                row++;
        }

        std::set<Point> rotatedCurves;
        Eigen::Affine3d rot(Eigen::AngleAxisd(angle, vectorRot));
        Eigen::Matrix<double, Eigen::Dynamic, 4> m = matrix * rot.matrix();
        for (int i = 0; i < m.innerSize(); i++) {
                if (m(i, 0) != 0 || m(i, 1) != 0 || m(i, 2) != 0) {
                        Scalar posx = m(i, 0);
                        Scalar posy = m(i, 1);
                        Scalar posz = m(i, 2) + 25;
                        rotatedCurves.insert({posx, posy, posz});
                }
        }

        std::set<Point> set;
        for (auto it = rotatedCurves.begin(), ite = rotatedCurves.end();
             it != ite; ++it) {
                Ball ball(*it, ballRadius);
                Point current = *it;
                std::vector<Point> pointsInBall = ball.pointsInBall();
                // drawEllipsoid(pointsInBall, 7, 10, 7, current[0], current[1], current[2], 0.5);
                set.insert(pointsInBall.begin(), pointsInBall.end());
        }
        Domain domain = PointUtil::computeBoundingBox<Domain>(set);
        Container aSet(domain);
        aSet.insert(set.begin(), set.end());
        return aSet;
}

template <typename Container>
Container Modeller<Container>::addNoise(const Container& set, float noise) {
        typedef Container Predicate;
        typedef DGtal::KanungoNoise<Predicate, Domain, Container> KanungoNoise;
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, Point::dimension> L2Metric;
        typedef DGtal::MetricAdjacency<Space, 3> MetricAdjacency;
        L2Metric l2Metric;
        KanungoNoise kanungo(set, set.domain(), noise);
        Container set2(set.domain());

        Container neighbors(set.domain());
        DGtal::DigitalSetInserter<Container> inserter(neighbors);
        for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
                MetricAdjacency::writeNeighbors(inserter, *it);
        }
        for (auto it = neighbors.begin(), itE = neighbors.end(); it != itE; ++it) {
                if (kanungo(*it)) {
                        set2.insertNew(*it);
                }
        }
        return set2;
}

template <typename Container>
Container Modeller<Container>::create2DCurve() {
        std::set<Point> digitalSet;
        digitalSet.insert({1,3});
        digitalSet.insert({2,3});
        digitalSet.insert({3,4});
        digitalSet.insert({4,4});
        digitalSet.insert({5,5});
        digitalSet.insert({6,5});
        digitalSet.insert({7,6});
        digitalSet.insert({8,6});
        digitalSet.insert({9,6});
        Domain domain = PointUtil::computeBoundingBox<Domain>(digitalSet);
        Container aSet(domain);
        aSet.insert(digitalSet.begin(), digitalSet.end());
        return aSet;
}

template <typename Container>
template <typename Board>
void Modeller<Container>::create2DNaiveTangentsForVisu(const Container & points, Board& board) {
        board << points.domain() << points;
        auto nextIt = points.begin();
        for (auto it = points.begin(),
                     itE = points.end();
             (++nextIt) != itE; ++it) {
                board.drawArrow((*it)[0], (*it)[1], (*nextIt)[0], (*nextIt)[1]);
        }
}




#endif
