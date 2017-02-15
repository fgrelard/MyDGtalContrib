#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <ShapeDescriptor.h>
#include <DGtal/io/writers/PGMWriter.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "geometry/DistanceToMeasure.h"
#include "shapes/Ball.h"
#include "geometry/SetProcessor.h"


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
        Z2i::RealVector v = delta.projection(p);
        Z2i::Point dest = p+v;
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
    for (const Z2i::RealVector& first : dirs) {
        int cptSame = 0, cptDifferent = 0;
        Z2i::RealVector f = first.getNormalized();
        for (const Z2i::RealVector& other : dirs) {
            Z2i::RealVector o = other.getNormalized();
            if (first == other) continue;
            if (f.dot(o) < 0)
                cptDifferent++;
            if (f.dot(o) >= 0)
                cptSame++;
        }
        double ratio = std::min(cptDifferent, cptSame) *1.0 / std::max(cptDifferent, cptSame);
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
        Z2i::RealVector vecProj = delta.projection(current);
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
        float v = ((float)*it) * 1.0 / thresholdMax;
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


    std::map<Z2i::Point, std::vector<Z2i::RealVector> > pToV = computePointToVectors(delta);
    std::map< Z2i::Point, double > pToValues = maxProjection(pToV, delta);
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();




/*
    DGtal::trace.beginBlock("Tube criterion");
    std::map<Z2i::Point, double> pToValues;
    for ( typename Domain::ConstIterator it = delta.domain().begin(),
                  itE = delta.domain().end(); it != itE; ++it ) {
        Point p = *it;
        Z3i::Point c(p[0], p[1], 0);
        float v = sqrt( delta.distance2( p ) );
        v = std::min( (float)m, std::max( v, 0.0f ) );
        int radius = 3;
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
//        if (angle > 0.52)
//            orientation = 0;
        double distanceNormalized = v / m;
        length = distanceNormalized;
//        double angle = twoAngleClusters(dirs);
//        length *= angle;
//        length = angle;
        pToValues[p] = length;
        RealVector grad = delta.projection(p);
        Z3i::RealVector grad3D(grad[0], grad[1], 0);
        if (std::isnan(grad3D.norm()) || grad3D.norm() == 0) continue;
        viewer.addLine(c, c+grad3D.getNormalized()*10);
    }
    DGtal::trace.endBlock()*/

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

    for ( const std::pair<Z2i::Point, double>& pair: pToValues) {
        Point p = pair.first;
        Z3i::Point c(p[0], p[1], 0);
        double value = pair.second;
        Color currentColor = cmap_grad(value);
        currentColor.alpha(180);
        RealVector grad = delta.projection(p);
        viewer << CustomColors3D(currentColor, currentColor)
               << c;

        Z3i::RealVector grad3D(grad[0], grad[1], 0);
        if (std::isnan(grad3D.norm()) || grad3D.norm() == 0) continue;
        viewer.addLine(c, c+grad3D);

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
