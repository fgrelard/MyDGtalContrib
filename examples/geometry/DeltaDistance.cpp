#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <ShapeDescriptor.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/geometry/volumes/distance/VoronoiMap.h>
#include <shapes/Polygon.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "geometry/DistanceToMeasure.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

class PointToVectors {
public:
    PointToVectors(const Z2i::Point& p) : myPoint(p) {}
    PointToVectors(const Z2i::Point& p,
                   const std::vector<Z2i::RealVector>& vectors) : myPoint(p), myVectors(vectors) {}

public:
    Z2i::Point point() const { return myPoint; }
    std::vector<Z2i::RealVector> vectors() { return myVectors; }
    void add(const Z2i::RealVector& vec) { myVectors.push_back(vec); }

    bool operator<(const PointToVectors &rhs) const {
        return myPoint < rhs.myPoint;
    }


private:
    Z2i::Point myPoint;
    std::vector<Z2i::RealVector> myVectors;
};

template <typename DistanceFunction>
std::map<Z2i::Point, std::vector<Z2i::RealVector> >
computePointToVectors(const DistanceFunction& delta) {
    std::map<Z2i::Point, std::vector<Z2i::RealVector> > aMap;
    for (const Z2i::Point& p : delta.domain()) {
        aMap[p] = std::vector<Z2i::RealVector>();
    }
    for (const Z2i::Point& p : delta.domain()) {
        Z2i::RealVector v = delta.projectionDistance(p);
        Z2i::Point dest = p+v;
        if (delta.domain().isInside(dest))
            aMap[dest].push_back(-v);
    }
    return aMap;
}

double signedAngle(const Z2i::RealVector& v) {
    double angle = atan2(v[1], v[0]);
    return (angle > 0 ? angle : (2*M_PI + angle));
}

double signedAngleTwoVectors(const Z2i::RealVector& v1,
                             const Z2i::RealVector& v2) {
    if (v1 == v2) return 2 * M_PI;
    double angle = signedAngle(v1) - signedAngle(v2);
    if (angle < 0) {
        return 2 * M_PI;
    }
    return angle;
}

double twoAngleClusters(const std::vector<Z2i::RealVector>& dirs) {
    if (dirs.size() == 1) return 0;
    double maxAngle = 0;
    for (const Z2i::RealVector& f : dirs) {
        Z2i::RealVector other = *min_element(dirs.begin(),
                                             dirs.end(),
                                             [&](const Z2i::RealVector& lhs,
                                                 const Z2i::RealVector& rhs) {
                                                 return (signedAngleTwoVectors(f, lhs) < signedAngleTwoVectors(f, rhs));
                                             });

        double angle = signedAngleTwoVectors(f, other);
        if (angle != 2* M_PI && angle > maxAngle)
            maxAngle = angle;
    }
    return maxAngle;
}

double meanAngle(const std::vector<Z2i::RealVector>& dirs) {
    std::vector<double> angles;
    for (const Z2i::RealVector& f : dirs) {
        Z2i::RealVector first = f.getNormalized();
        for (const Z2i::RealVector& o : dirs) {
            if (f==o) continue;
            Z2i::RealVector other = o.getNormalized();
            if (first.dot(other) < 0)
                other = -other;
            double angle = first.cosineSimilarity(other);
            angles.push_back(angle);
        }
    }
    if (angles.size() == 0) return 0.0;
    Statistic<double> statistic(true);
    statistic.addValues(angles.begin(), angles.end());
    double value = statistic.mean();
    return value;
}

double lengthVectors(const std::vector<Z2i::RealVector>& dirs) {
    if (dirs.size() == 1) return 0;
    Z2i::RealVector sum;
    for (const Z2i::RealVector& f : dirs) {
        sum += f;
    }
    return sum.norm() / dirs.size();
}

double twoOrientation(const std::vector<Z2i::RealVector>& dirs) {
    if (dirs.size() <= 2) return 0.0;
    std::vector<double> orientation;
    for (int i = 0; i < dirs.size(); i++) {
        Z2i::RealVector first = dirs[i];
        int cptSame = 0, cptDifferent = 0;
        Z2i::RealVector f = first.getNormalized();
        for (int j = 0; j < dirs.size(); j++) {
            if (i == j) continue;
            Z2i::RealVector other = dirs[j];
            Z2i::RealVector o = other.getNormalized();
            if (f.dot(o) < 0)
                cptDifferent++;
            else
                cptSame++;
        }
        int min = std::min(cptDifferent, cptSame);
        int max = std::max(cptDifferent, cptSame);
        double ratio = min *1.0 / max;
        orientation.push_back(ratio);
    }
    Statistic<double> stat(true);
    stat.addValues(orientation.begin(), orientation.end());
    double value = stat.median();
    double angle = meanAngle(dirs);
    value *= std::abs(M_PI/4 - angle);
    return value;
}

double lengthMinorAxis(const std::vector<Z2i::RealVector>& dirs) {
    Z2i::Domain domain(Z2i::Point(-10, -10), Z2i::Point(10,10));
    Z2i::DigitalSet aSet(domain);
    for (const Z2i::RealVector& v : dirs)
        aSet.insert(v*10);
    ShapeDescriptor<Z2i::DigitalSet> shapeDescriptor(aSet);
    auto matrix = shapeDescriptor.computeCovarianceMatrix();
    if (matrix.size() == 0)  return 0.0;
    double l0 = shapeDescriptor.extractEigenValue(matrix, 0);
    double l1 = shapeDescriptor.extractEigenValue(matrix, 1);
    return l1;
}

template <typename DistanceToMeasure>
std::map<Z2i::Point, double>
maxProjection(const std::map<Z2i::Point, std::vector<Z2i::RealVector> >& pToV,
              const DistanceToMeasure& delta) {
    std::map<Z2i::Point, double> aMap;
    for (const Z2i::Point& current:  delta.domain()) {
        std::vector<Z2i::RealVector> vec = pToV.at(current);
        aMap[current] = vec.size();
    }
    return aMap;
}

template <typename DistanceToMeasure>
std::map<Z2i::Point, double>
maxProjectionRadius(const std::map<Z2i::Point, std::vector<Z2i::RealVector> >& pToV,
                    const DistanceToMeasure& delta) {
    std::map<Z2i::Point, double> aMap;
    for (const Z2i::Point& current:  delta.domain()) {
        Z2i::RealVector vecProj = delta.projectionDistance(current);
        Z2i::Point projection = current  + vecProj;
        if (!delta.domain().isInside(projection)) continue;
        std::vector<Z2i::RealVector> vec = pToV.at(projection);
        Z2i::RealVector candidate = *std::max_element(vec.begin(), vec.end(), [&](const Z2i::RealVector& v1,
                                                                                  const Z2i::RealVector& v2) {
            return v1.norm() < v2.norm();
        });
        Z2i::Point proj = projection + candidate;
        double value = sqrt(delta.distance2(proj));
        aMap[current] = value;
    }
    return aMap;

}

template <typename Domain>
std::map<Z2i::Point, std::set<Z2i::Point> >
computeVoronoiMap(const std::map<Z2i::Point, std::vector<Z2i::RealVector> >& vectors,
                  const Domain& domainImage) {
    std::map<Z2i::Point, std::set<Z2i::Point> > aMap;
    for (const auto& pair : vectors) {
        Z2i::Point p = pair.first;
        std::vector<Z2i::RealVector> vec = pair.second;

        std::vector<Z2i::Point> vertices;

        for (const Z2i::RealVector& v : vec) {
            Z2i::Point proj = p + v;
            if (domainImage.isInside(proj))
                vertices.push_back(proj);
        }
        vertices.push_back(p);

        Z2i::Domain domain = PointUtil::computeBoundingBox<Z2i::Domain>(vertices);
        Z2i::Point lower = domain.lowerBound() - Z2i::Point::diagonal(1);
        Z2i::Point upper = domain.upperBound() + Z2i::Point::diagonal(1);
        domain = Domain(lower, upper);
        Polygon<Z2i::Point> polygon(vertices.begin(), vertices.end());
        std::set<Z2i::Point> hullSet;
        for (const Z2i::Point &d : domain) {
            if (polygon.isInside(d))
                hullSet.insert(d);
        }
        if (hullSet.size() != 0) {
            aMap[p] = hullSet;
        }
    }
    return aMap;
}

int main( int argc, char **argv )
{
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage2D;
    typedef ImageContainerBySTLVector<Domain, float>         FloatImage2D;
    typedef DistanceToMeasure<FloatImage2D>                 Distance;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
            ("mass,a", po::value<double>()->default_value(1), "Mass to integrate for distance to measure")
            ("rmax,r", po::value<double>()->default_value(10), "Max radius for delta distance")
            ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
            ;


    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, (const char *const *) argv, general_opt), vm);
    } catch(const std::exception& ex){
        parseOK=false;
        trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
    po::notify(vm);
    if( !parseOK || vm.count("help")||argc<=1)
    {
        std::cout << "Usage: " << argv[0] << " [input]\n"
                  << "Display volume file as a voxel set by using QGLviewer"<< endl
                  << general_opt << "\n";
        return 0;
    }
    if(!vm.count("input"))
    {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }


    string inputFilename = vm["input"].as<std::string>();
    int thresholdMax = vm["thresholdMax"].as<int>();
    double mass = vm["mass"].as<double>();
    double rmax = vm["rmax"].as<double>();

    GrayLevelImage2D img  = GenericReader<GrayLevelImage2D>::import( inputFilename );
    auto domain = img.domain();
    FloatImage2D     fimg( img.domain() );
    FloatImage2D::Iterator outIt = fimg.begin();
    for ( GrayLevelImage2D::ConstIterator it = img.begin(), itE = img.end(); it != itE; ++it ) {
        float v = ((float) *it) * 1.0f / thresholdMax;
        *outIt++ = v;
    }
    trace.beginBlock( "Computing delta-distance." );
    Distance     delta( mass, fimg, rmax );
    trace.endBlock();

    float m = 0.0f;
    for ( typename Domain::ConstIterator it = delta.domain().begin(),
                  itE = delta.domain().end(); it != itE; ++it )
    {
        Point p = *it;
        float v = sqrt( delta.distance2( p ) );
        m = std::max( v, m );
    }


    srand(time(NULL));
    std::map<Z2i::Point, std::vector<Z2i::RealVector> > pToV = computePointToVectors(delta);
    //std::map< Z2i::Point, double > pToValues = maxProjectionRadius(pToV, delta);
    //std::map< Z2i::Point, double > pToValues = maxProjection(pToV, delta);
    std::map<Z2i::Point, std::set<Z2i::Point> > voro = computeVoronoiMap(pToV, domain);
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();

    DGtal::trace.beginBlock("Tube criterion");
    std::map<Z2i::Point, double> pToValues;
    for ( typename Domain::ConstIterator it = delta.domain().begin(),
                  itE = delta.domain().end(); it != itE; ++it ) {
        Point p = *it;
        Z3i::Point c(p[0], p[1], 0);
        float v = sqrt( delta.distance2( p ) );
        v = std::min( (float)m, std::max( v, 0.0f ) );
        int radius = 4;
        Ball<Point> ball(p, radius);
        FloatImage2D aSet = ball.intersection( fimg );
        std::vector<RealVector> dirs;
        for (const Point& b : aSet.domain()) {
            if (aSet(b) == 0) continue;
            std::vector<RealVector> vec = pToV[b];
            for (const RealVector& v : vec) {
                if (//ball.contains(b+grad) &&
                        v != RealVector::zero) {
                    RealVector grad = v.getNormalized();
                    dirs.push_back(grad);
                }
            }

        }
        //No normalization
        double length = lengthMinorAxis(dirs);
        double angle = meanAngle(dirs);
        double orientation = twoOrientation(dirs);
        double distanceNormalized = v / m;
        length = v;
        pToValues[p] = length;

    }
    DGtal::trace.endBlock();

    double min = std::min_element(pToValues.begin(), pToValues.end(), [&](const std::pair<Z2i::Point, double>& p1,
                                                                          const std::pair<Z2i::Point, double>& p2) {
        return p1.second < p2.second;
    })->second;

    double max = std::max_element(pToValues.begin(), pToValues.end(), [&](const std::pair<Z2i::Point, double>& p1,
                                                                          const std::pair<Z2i::Point, double>& p2) {
        return p1.second < p2.second;
    })->second;
    DGtal::trace.info() << min << " " << max << std::endl;
    GradientColorMap<float> cmap_grad( min, max );

    cmap_grad.addColor( Color( 255, 255, 255 ) );
    cmap_grad.addColor( Color( 255, 255, 0 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 0, 255, 0 ) );
    cmap_grad.addColor( Color( 0,   0, 255 ) );
    cmap_grad.addColor( Color( 0,   0, 0 ) );

    for (const Z2i::Point &p : img.domain()) {
        Z3i::Point p3D(p[0], p[1], 0);
        if (img(p) == 255)
            viewer << CustomColors3D(Color::Red, Color::Red) << p3D;
    }


    for (const auto& pair : pToV)  {
        Point p = pair.first;
        std::vector<Z2i::RealVector> vec = pair.second;
        if (pToValues.find(p) == pToValues.end()) continue;
        Z3i::Point p3D(p[0], p[1], 0);
        double value = pToValues.at(p);
        Color currentColor = cmap_grad(value);
        currentColor.alpha(180);
        viewer << CustomColors3D(currentColor, currentColor) << p3D;

//        if (p[0] % 30 != 0 || p[1] % 30 != 0) continue;

        for (const Z2i::RealVector &v : vec) {
            Z3i::RealVector vec3D(v[0], v[1], 0);
            if (std::isnan(v.norm()) || v.norm() == 0) continue;
            viewer.addLine(p3D, p3D + vec3D);
        }
        if (voro.find(p) == voro.end()) continue;



        size_t size = vec.size();
        std::set<Z2i::Point> vorocell = voro.at(p);
        int r = rand() % 256, g = rand() % 256, b = rand() % 256;
        Color color(r, g, b);
        for (const Z2i::Point &v : vorocell) {
            Z3i::Point v3D(v[0], v[1], 0);
            viewer << CustomColors3D(color, color) << v3D;
        }

    }


    GrayLevelImage2D out(domain);
    for (const std::pair<Z2i::Point, double>& pair : pToValues) {
        unsigned char value = pair.second * 255 / max;
        out.setValue(pair.first, value);
    }

    PGMWriter<GrayLevelImage2D>::exportPGM("test.pgm", out);
    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
