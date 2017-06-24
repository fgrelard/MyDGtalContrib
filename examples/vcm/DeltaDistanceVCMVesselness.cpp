#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <ShapeDescriptor.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/geometry/volumes/distance/VoronoiMap.h>
#include <shapes/Polygon.h>
#include <geometry/DistanceToMeasure.h>
#include <viewer/ViewerDistanceBall.h>

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "geometry/DistanceToMeasureEdge.h"
#include "DGtal/io/readers/ITKReader.h"
#include "vcm/VCMAdjustableRadius.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "DGtal/io/writers/ITKWriter.h"
#include "vcm/VCMVesselness.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

class PointToVectors {
public:
    PointToVectors(const Z3i::Point& p) : myPoint(p) {}
    PointToVectors(const Z3i::Point& p,
                   const std::vector<Z3i::RealVector>& vectors) : myPoint(p), myVectors(vectors) {}

public:
    Z3i::Point point() const { return myPoint; }
    std::vector<Z3i::RealVector> vectors() { return myVectors; }
    void add(const Z3i::RealVector& vec) { myVectors.push_back(vec); }

    bool operator<(const PointToVectors &rhs) const {
        return myPoint < rhs.myPoint;
    }


private:
    Z3i::Point myPoint;
    std::vector<Z3i::RealVector> myVectors;
};

template <typename DistanceFunction>
std::map<Z3i::Point, std::vector<Z3i::RealVector> >
computePointToVectors(const DistanceFunction& delta) {
    std::map<Z3i::Point, std::vector<Z3i::RealVector> > aMap;
    for (const Z3i::Point& p : delta.domain()) {
        aMap[p] = std::vector<Z3i::RealVector>();
    }
    for (const Z3i::Point& p : delta.domain()) {
        Z3i::RealVector v = delta.projectionDistance(p);
        Z3i::Point dest = p+v;
        if (delta.domain().isInside(dest))
            aMap[dest].push_back(-v);
    }
    return aMap;
}

template <typename DistanceFunction>
std::map<Z3i::Point, std::vector<Z3i::RealVector> >
computePointToVectorsSobel(const DistanceFunction& delta) {
    std::map<Z3i::Point, std::vector<Z3i::RealVector> > aMap;
    for (const Z3i::Point& p : delta.domain()) {
        aMap[p] = std::vector<Z3i::RealVector>();
    }
    for (const Z3i::Point& p : delta.domain()) {
        Z3i::RealVector v = delta.projection(p);
        if (v == Z3i::RealVector::zero) {
            aMap[p].push_back(Z3i::RealVector::zero);
            continue;
        }
        v = v.getNormalized();
        Z3i::RealPoint dest = (Z3i::RealPoint)p;
        Z3i::RealPoint next = dest + v;
        while (delta.domain().isInside(next) && delta((Z3i::Point)next) <= delta((Z3i::Point)dest)) {
            dest = next;
            next = dest + v;
        }
        Z3i::Point candidate = dest;
        if (!delta.domain().isInside(next)) {
            aMap[candidate].push_back(Z3i::RealVector::zero);
        }
        else if (delta.domain().isInside(candidate))
            aMap[candidate].push_back(-v);
    }
    return aMap;
}

double signedAngle(const Z3i::RealVector& v) {
    double angle = atan2(v[1], v[0]);
    return (angle > 0 ? angle : (2*M_PI + angle));
}

double signedAngleTwoVectors(const Z3i::RealVector& v1,
                             const Z3i::RealVector& v2) {
    if (v1 == v2) return 2 * M_PI;
    double angle = signedAngle(v1) - signedAngle(v2);
    if (angle < 0) {
        return 2 * M_PI;
    }
    return angle;
}

double twoAngleClusters(const std::vector<Z3i::RealVector>& dirs) {
    if (dirs.size() == 1) return 0;
    double maxAngle = 0;
    for (const Z3i::RealVector& f : dirs) {
        Z3i::RealVector other = *min_element(dirs.begin(),
                                             dirs.end(),
                                             [&](const Z3i::RealVector& lhs,
                                                 const Z3i::RealVector& rhs) {
                                                 return (signedAngleTwoVectors(f, lhs) < signedAngleTwoVectors(f, rhs));
                                             });

        double angle = signedAngleTwoVectors(f, other);
        if (angle != 2* M_PI && angle > maxAngle)
            maxAngle = angle;
    }
    return maxAngle;
}

double meanAngle(const std::vector<Z3i::RealVector>& dirs) {
    std::vector<double> angles;
    for (const Z3i::RealVector& f : dirs) {
        Z3i::RealVector first = f.getNormalized();
        for (const Z3i::RealVector& o : dirs) {
            if (f==o) continue;
            Z3i::RealVector other = o.getNormalized();
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

double lengthVectors(const std::vector<Z3i::RealVector>& dirs) {
    if (dirs.size() == 1) return 0;
    Z3i::RealVector sum;
    for (const Z3i::RealVector& f : dirs) {
        sum += f.getNormalized();
    }
    return sum.norm();
}

double twoOrientation(const std::vector<Z3i::RealVector>& dirs) {
    if (dirs.size() <= 2) return 0.0;
    std::vector<double> orientation;
    for (int i = 0; i < dirs.size(); i++) {
        Z3i::RealVector first = dirs[i];
        int cptSame = 0, cptDifferent = 0;
        Z3i::RealVector f = first.getNormalized();
        for (int j = 0; j < dirs.size(); j++) {
            if (i == j) continue;
            Z3i::RealVector other = dirs[j];
            Z3i::RealVector o = other.getNormalized();
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
    //value *= std::abs(M_PI/4 - angle);
    return value;
}

double lengthMinorAxis(const std::vector<Z3i::RealVector>& dirs) {
    Z3i::Domain domain(Z3i::Point(-10, -10), Z3i::Point(10,10));
    Z3i::DigitalSet aSet(domain);
    for (const Z3i::RealVector& v : dirs)
        aSet.insert(v*10);
    ShapeDescriptor<Z3i::DigitalSet> shapeDescriptor(aSet);
    auto matrix = shapeDescriptor.computeCovarianceMatrix();
    if (matrix.size() == 0)  return 0.0;
    double l0 = shapeDescriptor.extractEigenValue(matrix, 0);
    double l1 = shapeDescriptor.extractEigenValue(matrix, 1);
    return l1;
}

bool areVectorsValid(const std::vector<Z3i::RealVector>& dirs) {
    for (const Z3i::RealVector& v1 : dirs) {
        Z3i::RealVector v1norm = v1.getNormalized();
        for (const Z3i::RealVector& v2 : dirs) {
            Z3i::RealVector v2norm = v2.getNormalized();
            if (v1norm.dot(v2norm) < 0)
                return false;
        }
    }
    return true;
}

template <typename DistanceToMeasure>
std::map<Z3i::Point, double>
maxProjection(const std::map<Z3i::Point, std::vector<Z3i::RealVector> >& pToV,
              const DistanceToMeasure& delta) {
    std::map<Z3i::Point, double> aMap;
    for (const Z3i::Point& current:  delta.domain()) {
        std::vector<Z3i::RealVector> vec = pToV.at(current);
        aMap[current] = vec.size();
    }
    return aMap;
}

template <typename DistanceToMeasure>
std::map<Z3i::Point, double>
maxProjectionRadius(const std::map<Z3i::Point, std::vector<Z3i::RealVector> >& pToV,
                    const DistanceToMeasure& delta) {
    std::map<Z3i::Point, double> aMap;
    for (const Z3i::Point& current:  delta.domain()) {
        Z3i::RealVector vecProj = delta.projectionDistance(current);
        Z3i::Point projection = current  + vecProj;
        if (!delta.domain().isInside(projection)) continue;
        std::vector<Z3i::RealVector> vec = pToV.at(projection);
        Z3i::RealVector candidate = *std::max_element(vec.begin(), vec.end(), [&](const Z3i::RealVector& v1,
                                                                                  const Z3i::RealVector& v2) {
                                                          return v1.norm() < v2.norm();
                                                      });
        double value = candidate.norm();
        aMap[current] = value;
    }
    return aMap;

}

template <typename Domain>
std::map<Z3i::Point, std::set<Z3i::Point> >
computeVoronoiMap(const std::map<Z3i::Point, std::vector<Z3i::RealVector> >& vectors,
                  const Domain& domainImage) {
    std::map<Z3i::Point, std::set<Z3i::Point> > aMap;
    for (const auto& pair : vectors) {
        Z3i::Point p = pair.first;
        std::vector<Z3i::RealVector> vec = pair.second;

        std::vector<Z3i::Point> vertices;

        for (const Z3i::RealVector& v : vec) {
            Z3i::Point proj = p + v;
            if (domainImage.isInside(proj))
                vertices.push_back(proj);
        }
        vertices.push_back(p);

        Z3i::Domain domain = PointUtil::computeBoundingBox<Z3i::Domain>(vertices);
        Z3i::Point lower = domain.lowerBound() - Z3i::Point::diagonal(1);
        Z3i::Point upper = domain.upperBound() + Z3i::Point::diagonal(1);
        domain = Domain(lower, upper);
        Polygon<Z3i::Point> polygon(vertices.begin(), vertices.end());
        std::set<Z3i::Point> hullSet;
        for (const Z3i::Point &d : domain) {
            if (polygon.isInside(d))
                hullSet.insert(d);
        }
        if (hullSet.size() != 0) {
            aMap[p] = hullSet;
        }
    }
    return aMap;
}

template <typename FImg>
FImg invert(const FImg &img) {
    FImg invert(img.domain());
    for (const auto &p : img.domain()) {
        invert.setValue(p, 255 - img(p));
    }
    return invert;
}

template <typename Matrix>
Matrix sqrtMatrix(const Matrix& in) {
    Matrix out;
    for (size_t i = 0; i < in.rows(); i++) {
        for (size_t j = 0; j < in.cols(); j++) {
            auto valueIJ = in(i,j);
            if (std::isnan(valueIJ)) out(i,j) = 0;
            if (valueIJ >= 0) {
                out(i, j) = std::sqrt(valueIJ);
            } else {
                out(i, j) = - std::sqrt(std::abs(valueIJ));
            }
        }
    }
    return out;
}


void exportToSDP(const
                 std::map<DGtal::Z3i::Point,
                 std::vector<DGtal::Z3i::RealVector> >& pointToVectors,
                 const std::string& outputFilename) {
    std::ofstream outStream;
    outStream.open(outputFilename.c_str());
    for (const auto & pToV : pointToVectors) {
        DGtal::Z3i::Point p = pToV.first;
        std::vector<DGtal::Z3i::RealVector> vectors = pToV.second;
        for (const DGtal::Z3i::RealVector& v : vectors) {
            if (v == DGtal::Z3i::RealVector::zero) continue;
            DGtal::Z3i::RealPoint dest = (DGtal::Z3i::RealPoint)p - v;
            DGtal::Z3i::RealVector vNorm = -v.getNormalized();
            outStream << p[0] << " " << p[1] << " " << p[2] << " " << dest[0] << " " << dest[1] << " " << dest[2] << " " << vNorm[0] << " " << vNorm[1] << " " <<vNorm[2] << std::endl;
        }
    }
    outStream.close();
}

int main( int argc, char **argv )
{
    using namespace DGtal;
    using namespace DGtal::Z3i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage;
    typedef ImageContainerBySTLVector<Domain, float>         FloatImage;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef DistanceToMeasureEdge<FloatImage> Distance;
    typedef VCMAdjustableRadius<Space, L2Metric> VCM;
    typedef functors::BallConstantPointFunction<Point, double> KernelFunction;
    typedef EigenDecomposition<3, double> LinearAlgebraTool;
    typedef PointVector<3, double> RealVector2f;
    typedef typename VCM::MatrixNN Matrix;
    typedef ImageContainerBySTLVector<Domain, Matrix> MatrixImage;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("distance,d", po::value<std::string>(), "vol file (corresponding volume)")
        ("scale,s", po::value<std::string>(), "scale image for VCM radius")
        ("output,o", po::value<std::string>(), "vol file (corresponding volume)")
        ("alpha,a", po::value<double>()->default_value(0.1), "Frangi alpha")
        ("beta,b", po::value<double>()->default_value(1), "Frangi beta")
        ("gamma,c", po::value<double>()->default_value(10), "Frangi gamma")
        ("radiusInside,n", po::value<double>()->default_value(10))
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
    if(!vm.count("distance"))
    {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }



    string distanceFilename = vm["distance"].as<std::string>();
    string scaleImage;
    if (vm.count("scale"))
        scaleImage = vm["scale"].as<std::string>();
    string outputFilename = vm["output"].as<std::string>();
    int thresholdMax = vm["thresholdMax"].as<int>();
    double alpha = vm["alpha"].as<double>();
    double beta = vm["beta"].as<double>();
    double gamma = vm["gamma"].as<double>();
    double radiusVCM = vm["radiusInside"].as<double>();


    FloatImage fimg = ITKReader<FloatImage>::importITK( distanceFilename );
    FloatImage fimgScale(fimg.domain());
    if (vm.count("scale"))
        fimgScale =  ITKReader<FloatImage>::importITK(scaleImage);

    auto domain = fimg.domain();


    trace.beginBlock( "Computing delta-distance." );
    Distance     delta ( fimg );
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
    std::map<Z3i::Point, std::vector<Z3i::RealVector> > pToV = computePointToVectors(delta);
    string newfilename= outputFilename.substr(0,outputFilename.find_last_of('.'))+".sdp";
    exportToSDP(pToV, newfilename);

    MatrixImage* imageMatrix = new MatrixImage(fimg.domain());
    DGtal::trace.beginBlock("Constructing delta VCM");
    for ( typename Domain::ConstIterator it = delta.domain().begin(),
              itE = delta.domain().end(); it != itE; ++it ) {
        Point p = *it;
        if (fimg(p) != 0)  {
            Matrix mat;
            imageMatrix->setValue(p, mat);
            continue;
        }
        if (vm.count("scale") && fimgScale(p) > 0)
            radiusVCM = fimgScale(p) + sqrt(3);

        Ball<Point> ball(p, radiusVCM);
        FloatImage aSet = ball.intersection( fimg );
        std::vector<RealVector> dirs;
        double sumDistanceB = 0, sumDistanceBV = 0;
        bool reverse = false;
        for (const Point& b : aSet.domain()) {
            if (aSet(b) == 0) continue;
            std::vector<RealVector> vec = pToV.at(b);
            for (const RealVector& v : vec) {
                if (v != RealVector::zero) {
                    RealVector grad = v;
                    RealPoint d = b + grad;
                    RealVector bp = b - p;
                    RealVector dp = d - p;
                    double distanceB  = Z3i::l2Metric(p, b);
                    double distanceBV = Z3i::l2Metric(p, b+v);
                    sumDistanceB += distanceB;
                    sumDistanceBV += distanceBV;
                    grad *= distanceB;
                    if ( bp.dot(dp) > 0)
                        dirs.push_back(grad);
                    else {
                        reverse = true;
                        break;
                    }
                }
            }
        }

        if (dirs.size() == 0 // || sumDistanceB > sumDistanceBV || reverse
            ) {
            Matrix mat;
            imageMatrix->setValue(p, mat);
        }
        else {
            Matrix mat;
            for (const RealVector& v : dirs) {
                Matrix tmp;
                for ( Dimension  i = 0; i < 3; ++i )
                    for ( Dimension j = 0; j < 3; ++j )
                        tmp.setComponent( i, j, v[ i ] * v[ j ] );
                mat += tmp;
            }
            mat = sqrtMatrix(mat);
            imageMatrix->setValue(p, mat);
        }
    }
    DGtal::trace.endBlock();


    DGtal::trace.beginBlock("Vesselness");
    VCMVesselness<MatrixImage> vesselness(shared_ptr<MatrixImage>(imageMatrix), alpha, beta, gamma);
    auto vesselnessImage = vesselness.computeVesselness();
    DGtal::trace.endBlock();

    FloatImage out(vesselnessImage.domain());
    for (const Point& p : out.domain()) {
        double value = vesselnessImage(p);
        if (std::isnan(value))
            out.setValue(p, 0);
        else {
            out.setValue(p, value);
        }
    }
    ITKWriter<FloatImage>::exportITK(outputFilename, out);
    return 0;

}
