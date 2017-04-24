//
// Created by florent on 13/04/17.
//

#include <iostream>
#include <fstream>
#include <sstream>
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

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    using namespace DGtal;
    using namespace DGtal::Z3i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> Image;
    typedef ImageContainerBySTLVector<Domain, float> FloatImage;
    typedef VCMAdjustableRadius<Space, L2Metric> VCM;
    typedef functors::BallConstantPointFunction<Point, double> KernelFunction;
    typedef EigenDecomposition<3, double> LinearAlgebraTool;
    typedef PointVector<3, double> RealVector2f;
    typedef typename VCM::MatrixNN Matrix;
    typedef ImageContainerBySTLVector<Domain, Matrix> MatrixImage;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
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
    int thresholdMax = vm["thresholdMax"].as<int>();
    int thresholdMin = vm["thresholdMin"].as<int>();
    double alpha = vm["alpha"].as<double>();
    double beta = vm["beta"].as<double>();
    double gamma = vm["gamma"].as<double>();
    double R = vm["radiusInside"].as<double>();
    double r = vm["radiusNeighbour"].as<double>();

    Image volume = ITKReader<Image>::importITK(inputFilename);
    Z3i::Domain domainVolume = volume.domain();
    Z3i::DigitalSet setVolume(domainVolume);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                  thresholdMin-1, thresholdMax);

    KernelFunction chi(1.0, r);
    L2Metric l2Metric;
    VCM vcm(R, r, l2Metric, true);
    vcm.init(setVolume.begin(), setVolume.end());
    Matrix vcm_r, evec, nil;
    RealVector2f eval;
    MatrixImage matrixImage(domainVolume);
    DGtal::trace.beginBlock("Matrix image");
    for (const Point& p : matrixImage.domain()) {
        matrixImage.setValue(p, nil);
    }
    for (const Point& p : setVolume) {
        vcm_r = vcm.measure(chi, p);
        if (vcm_r == nil) continue;
        Matrix tmp;
        for (size_t i = 0; i < vcm_r.rows(); i++) {
            for (size_t j = 0; j < vcm_r.cols(); j++) {
                auto valueIJ = vcm_r(i,j);
                if (valueIJ >= 0) {
                    tmp(i, j) = std::sqrt(valueIJ);
                } else {
                    tmp(i, j) = - std::sqrt(std::abs(valueIJ));
                }
            }
        }
        matrixImage.setValue(p, tmp);
    }
    DGtal::trace.endBlock();

    DGtal::trace.beginBlock("Vesselness");
    VCMVesselness<MatrixImage> vesselness(matrixImage, alpha, beta, gamma);
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
    ITKWriter<FloatImage>::exportITK(outname, out);
    return 0;
}
