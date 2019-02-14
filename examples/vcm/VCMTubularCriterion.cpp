//
// Created by florent on 13/04/17.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/io/writers/ITKWriter.h>
#include <DGtal/io/readers/ITKReader.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include <algorithm>
#include <vcm/VCMAdjustableRadius.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <vcm/VCMVesselness.h>
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/math/Statistic.h"
#include "shapes/Ball.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "ShapeDescriptor.h"
#include "shapes/DigitalPlane.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

namespace DGtal {
    typedef SimpleMatrix<double, 2, 2> MatrixDouble;

    bool operator!=(const MatrixDouble &m1, const MatrixDouble &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<float, 2, 2> MatrixFloat;

    bool operator!=(const MatrixFloat &m1, const MatrixFloat &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<double, 3, 3> MatrixDouble3D;

    bool operator!=(const MatrixDouble3D &m1, const MatrixDouble3D &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<float, 3, 3> MatrixFloat3D;

    bool operator!=(const MatrixFloat3D &m1, const MatrixFloat3D &m2) { return !(m1 == m2); }

    namespace functors {
        bool operator==(Identity f1, Identity f2) { return true; }
    }
}

template <typename RealVector>
double eigenNorm(const RealVector& eval) {
    double l0 = std::sqrt(eval[0]);
    double l1 = std::sqrt(eval[1]);
    //double l2 = std::sqrt(eval[2]);
    double eigenNorm = std::sqrt(l0 * l0 + l1 * l1);
    return l1;
}


template <typename Plane, typename DigitalSet, typename Image, typename Point>
Plane planeForTube(const Image& image,
                 const std::map<Point, std::vector<Point> >& sToM,
                 const Point& site) {
    typedef typename DigitalSet::Space Space;
    typedef typename DigitalSet::Domain Domain;
    typedef DGtal::SimpleMatrix<double, Space::dimension, Space::dimension> Matrix;
    typedef EigenDecomposition<Space::dimension, double, Matrix> LinearAlgebraTool;
    typedef PointVector<Space::dimension, double> RealVector;

    std::vector<Point> cell = sToM.at(site);
    auto domain = PointUtil::computeBoundingBox<Domain>(cell);
    DigitalSet cellSet(image.domain());
    cellSet.insert(cell.begin(), cell.end());

    ShapeDescriptor<DigitalSet> sc(cellSet);
    auto matP = sc.computeCovarianceMatrix();
    if (matP.rows() == 0 || matP.cols() == 0) return Plane(Point::zero, RealVector::zero);
    Matrix dgMatP;
    for (int i = 0; i < Space::dimension; i++) {
        for (int j = 0; j < Space::dimension; j++) {
            dgMatP.setComponent(i, j, matP(i,j));
        }
    }

    Matrix evec;
    RealVector eval;
    LinearAlgebraTool::getEigenDecomposition(dgMatP, evec, eval);

    RealVector v2 = evec.column(Space::dimension-1);
    Plane plane(site, v2);

    DigitalSet above(image.domain()), below(image.domain());

    DigitalSet localArea(image.domain());
    Ball<Point> ball(site, 5);
    for (const Point& b : ball.pointSet()) {
        if (image.domain().isInside(b))
            localArea.insert(b);
    }
    for (const Point& p : localArea) {
        if (plane.isPointAbove(p))
            above.insert(p);
        else
            below.insert(p);
    }

    Statistic<double> statA(true), statB(true);

    for (const Point& p : above)
        statA.addValue(image(p));
    for (const Point& p : below)
        statB.addValue(image(p));

    if (statA.mean() > statB.mean())
        return plane;
    plane = Plane(site, -v2);
    return plane;
}

template <typename Point, typename DigitalSet>
bool isOutOfTube(const std::map<Point, std::vector<Point> >& sToM,
                 const Point& currentPoint,
                 const Point& site,
                 double radius) {
    typedef typename DigitalSet::Space Space;
    typedef typename DigitalSet::Domain Domain;
    typedef DGtal::SimpleMatrix<double, Space::dimension, Space::dimension> Matrix;
    typedef EigenDecomposition<Space::dimension, double, Matrix> LinearAlgebraTool;
    typedef PointVector<Space::dimension, double> RealVector;

    Statistic<double> X(true);


    Ball<Point> ballP(currentPoint, radius);
    Ball<Point> ballS(site, radius);

    std::vector<Point> cell = sToM.at(site);
    auto domain = PointUtil::computeBoundingBox<Domain>(cell);
    DigitalSet cellSet(domain);
    cellSet.insert(cell.begin(), cell.end());

    DigitalSet fromSites = ballS.intersection(cellSet);
    DigitalSet fromP = ballP.intersection(cellSet);

    // ShapeDescriptor<DigitalSet> sdp(fromP);
    // ShapeDescriptor<DigitalSet> sds(fromSites);

    // auto matP = sdp.computeCovarianceMatrix();
    // auto matS = sds.computeCovarianceMatrix();

    // Matrix dgMatP, dgMatS;

    // for (int i = 0; i < Space::dimension; i++) {
    //     for (int j = 0; j < Space::dimension; j++) {
    //         dgMatP.setComponent(i, j, matP(i,j));
    //         dgMatS.setComponent(i, j, matS(i,j));
    //     }
    // }

    // Matrix evec, evecS;
    // RealVector eval, evalS;
    // LinearAlgebraTool::getEigenDecomposition(dgMatP, evec, eval);
    // LinearAlgebraTool::getEigenDecomposition(dgMatS, evecS, evalS);

    // double nP = eigenNorm(eval);
    // double nS = eigenNorm(evalS);

    return (fromSites.size() < fromP.size());
}


template <typename DigitalSet, typename Point>
DigitalSet extractSites(const Point& p, double radius, const DigitalSet& sites) {
    Ball<Point> ball(p, radius);
    DigitalSet intersection = ball.intersection(sites);
    return intersection;
}

template <typename TImage>
TImage compute_vcm(const TImage& image, const TImage& grayImage,
                 std::string outname,
                 int thresholdMax, double alpha, double beta, double gamma, double R, double r) {
    typedef typename TImage::Domain Domain;
    typedef typename Domain::Space Space;
    typedef typename Space::Point Point;
    typedef ImageContainerBySTLVector<Domain, double> DoubleImage;
    typedef ImageContainerBySTLVector<Domain, float> FloatImage;
    typedef functors::HatPointFunction<Point, float> KernelFunction;
    typedef EigenDecomposition<Space::dimension, double> LinearAlgebraTool;
    typedef ImageContainerBySTLVector<Domain, unsigned char> Image;
    typedef ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef VCMAdjustableRadius<Space, L2Metric> VCM;
    typedef typename VCM::MatrixNN Matrix;
    typedef PointVector<Space::dimension, double> RealVector2f;
    typedef typename DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet;
    typedef functors::NotPointPredicate<DigitalSet> NotPredicate;
    typedef VoronoiMap<Space, NotPredicate, L2Metric > Voronoi;
    typedef DigitalPlane<Space> Plane;


//   Image volume = VolReader<Image>::importVol(inputFilename);
    Domain domainVolume = image.domain();
    DigitalSet setVolume(domainVolume);
    SetFromImage<DigitalSet>::append(setVolume, image,
                                     1, thresholdMax);
    DGtal::trace.info() << setVolume.size() << std::endl;
    L2Metric l2Metric;
    NotPredicate notSetPred(setVolume);
    Voronoi voronoimap(domainVolume,notSetPred,l2Metric);

    KernelFunction chi(1.0, r);
    VCM vcm(R, r, l2Metric, true);
    vcm.init(setVolume.begin(), setVolume.end());
    Matrix vcm_r, evec, nil;
    RealVector2f eval;
    FloatImage out(domainVolume);
    Image grayOut(domainVolume);
    DGtal::trace.beginBlock("Matrix image");

    int i = 0;


    std::map<Point, std::vector<Point> > siteToVCM;
    //Compute VCM for sites
    for (const Point& p : setVolume) {
        siteToVCM[p] = std::vector<Point>();

    }

     for (const Point& p : image.domain()) {
         Point site =voronoimap(p);
         siteToVCM[site].push_back(p);
     }


    std::map<Point, Plane> planes;
    for (const Point& p : setVolume) {
        Plane plane = planeForTube<Plane, DigitalSet>(grayImage, siteToVCM, p);
        planes[p] = plane;
    }

    vcm.setMySmallR(r);
    for (const Point& p : domainVolume) {
        DGtal::trace.progressBar(i, domainVolume.size());
        i++;
        // vcm_r = vcm.measure(chi, p);
        std::vector<Matrix> matrices = vcm.measureMedian(p);
        vcm_r = nil;
        for (Matrix m : matrices) {
            vcm_r += m;
        }
        if (vcm_r == nil) continue;


        Point site = voronoimap(p);

        //        DigitalSet sites = extractSites(p, r, setVolume);
        // double meanNormSites = meanEigenNorm(siteToVCM,sites);
        // double normP = eigenNorm(eval);

        // DGtal::trace.info() << meanNormSites << " " << normP << std::endl;

        // if (normP > meanNormSites) continue;
//        bool outOfTube = isOutOfTube<Point, DigitalSet>(siteToVCM, p, site, r);
//        if (outOfTube) continue;
        Plane plane = planes[site];
        if (plane.getCenter() != Point::zero && !plane.isPointAbove(p)) continue;


        LinearAlgebraTool::getEigenDecomposition(vcm_r, evec, eval);
        double l0 = std::sqrt(eval[0]);
        double l1 = std::sqrt(eval[1]);

        if (Space::dimension==2) {
            double circularity = l0 / l1;
            double noise = std::sqrt( l0 * l0 + l1 * l1);
            double criterion =
                exp(-(std::pow(circularity, 2) / (2 * std::pow(alpha, 2)))) *
                (1.0 - exp(-(std::pow(noise,2) / (2 * std::pow(gamma, 2)))));
            out.setValue(p, criterion);
        }
        if (Space::dimension == 3) {
            double norm = 1.0;
            if (matrices.size() > 0) {
                Matrix median = matrices[matrices.size()/2];
                norm = std::sqrt(median(0,0) + median(1,1) + median(2,2));
                norm = (norm == 0) ? 1.0 : norm;
            }
            double l2 = std::sqrt(eval[2]);
            double tub = l0 / std::sqrt(l1 * l2);
            double circularity = (l1 / l2) // * (1.0 / norm)
                ;

            double noise = std::sqrt( l0 * l0 + l1 * l1 + l2 * l2);
            double criterion = (1.0 - exp(-(std::pow(circularity, 2) / (2 * std::pow(alpha, 2))))) *
                exp(-(std::pow(tub, 2) / (2 * std::pow(beta, 2)))) *
                (1.0 - exp(-(std::pow(noise,2) / (2 * std::pow(gamma, 2)))));
            criterion = circularity;
            out.setValue(p, criterion);
            grayOut.setValue(p, (int)(criterion*255));

        }
    }
    DGtal::trace.endBlock();


    ITKWriter<FloatImage>::exportITK(outname, out);
    return grayOut;
}

int main(int argc, char **argv) {
    using namespace DGtal;



    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
        ("grayimage,g", po::value<std::string>(), "gray image (corresponding volume)")
        ("output,o", po::value<std::string>(), "vol file (corresponding volume)")
        ("alpha,a", po::value<double>()->default_value(0.1), "Frangi alpha")
        ("beta,b", po::value<double>()->default_value(1), "Frangi beta")
        ("gamma,c", po::value<double>()->default_value(10), "Frangi gamma")
        ("radiusInside,R", po::value<double>()->default_value(10), "radius of the ball inside voronoi cell")
        ("radiusNeighbour,n", po::value<double>()->default_value(10), "radius of the ball for the neighbourhood")
        ("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
        ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization");

    bool parseOK = true;
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    } catch (const std::exception &ex) {
        parseOK = false;
        trace.info() << "Error checking program options: " << ex.what() << endl;
    }
    po::notify(vm);
    if (!parseOK || vm.count("help") || argc <= 1) {
        std::cout << "Usage: " << argv[0] << " [input]\n"
                  << "Display volume file as a voxel set by using QGLviewer" << endl
                  << general_opt << "\n";
        return 0;
    }
    if (!vm.count("input")) {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }

    string inputFilename = vm["input"].as<std::string>();
    string outname = vm["output"].as<std::string>();
    string grayname = vm["grayimage"].as<std::string>();
    int thresholdMax = vm["thresholdMax"].as<int>();
    int thresholdMin = vm["thresholdMin"].as<int>();
    double alpha = vm["alpha"].as<double>();
    double beta = vm["beta"].as<double>();
    double gamma = vm["gamma"].as<double>();
    double R = vm["radiusInside"].as<double>();
    double r = vm["radiusNeighbour"].as<double>();


     const std::string extension =
      inputFilename.substr( inputFilename.find_last_of( "." ) + 1 );
     if ( extension == "pgm" || extension == "png" || extension == "gif" )
     { // go for 2D
         using namespace DGtal::Z2i;
         typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage;
          GrayLevelImage img = GenericReader<GrayLevelImage>::import(inputFilename);
          GrayLevelImage grayImg = GenericReader<GrayLevelImage>::import(grayname);
          GrayLevelImage grayOut = compute_vcm(img, grayImg, outname, thresholdMax,  alpha, beta, gamma, R, r);
     }
     else
     { // go for 3D
         using namespace DGtal::Z3i;
         typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage;
          GrayLevelImage img = GenericReader<GrayLevelImage>::import(inputFilename);
          GrayLevelImage grayImg = GenericReader<GrayLevelImage>::import(grayname);
          GrayLevelImage grayOut = compute_vcm(img, grayImg, outname, thresholdMax,  alpha, beta, gamma, R, r);
          std::string basename = boost::filesystem::basename(outname);
          VolWriter<GrayLevelImage>::exportVol(basename+".vol", grayOut);
     }

    return 0;
}
