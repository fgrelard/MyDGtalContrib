#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <ShapeDescriptor.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/geometry/volumes/distance/VoronoiMap.h>
#include <geometry/DistanceToMeasure.h>
#include <viewer/ViewerDistanceBall.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkDiscreteHessianGaussianImageFunction.h>
#include <DGtal/io/writers/ITKWriter.h>
#include <itkHessianToObjectnessMeasureImageFilter.h>
#include <itkMultiScaleHessianBasedMeasureImageFilter.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "geometry/DistanceToMeasureEdge.h"
#include "hessian/HessianRecursiveGaussian.h"
#include "hessian/FrangiVesselness.h"

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

template <typename FImg>
FImg invert(const FImg &img) {
    FImg invert(img.domain());
    for (const auto &p : img.domain()) {
        invert.setValue(p, 255 - img(p));
    }
    return invert;
}



template <typename TImage>
struct ImageToDistance {
        typedef typename TImage::Domain Domain;
        typedef typename Domain::Point Point;

        ImageToDistance() : myImage(Domain(Point::zero, Point::zero)), myDistance(std::numeric_limits<float>::max()) {}

    ImageToDistance(TImage image, float distance) : myImage(image), myDistance(distance) {}

    ImageToDistance(const ImageToDistance& other) : myImage(other.myImage), myDistance(other.myDistance) {}

    float myDistance;
    TImage myImage;
};

struct ValueToSqr {
    float operator()(float value) {
        return value * value;
    }
};

template <typename TImage>
TImage findClosestImage(const std::vector<ImageToDistance<TImage> >& vImage, float distance) {
    typedef ImageToDistance<TImage> Value;

    Value minImage;
    for (const Value& imageToDistance : vImage) {
        float diff = std::abs(imageToDistance.myDistance - distance);
        if (diff < minImage.myDistance) {
            minImage = imageToDistance;
        }
    }
    return minImage.myImage;
}

int main(int argc, char **argv) {
    using namespace DGtal;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "image file (corresponding volume)")
            ("distance,d", po::value<std::string>(), "distance image")
            ("output,o", po::value<std::string>(), "image file (corresponding volume)")
            ("mass,a", po::value<double>()->default_value(1), "Mass to integrate for distance to measure")
            ("rmax,r", po::value<double>()->default_value(10), "Max radius for delta distance")
            ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization");


    bool parseOK = true;
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, (const char *const *) argv, general_opt), vm);
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
    string distanceName;

    int thresholdMax = vm["thresholdMax"].as<int>();
    double mass = vm["mass"].as<double>();
    double rmax = vm["rmax"].as<double>();
    constexpr unsigned int dimension = 3;

    typedef SpaceND<dimension> Space;
    typedef HyperRectDomain<Space> Domain;
    typedef typename Space::Point Point;
    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage;

    typedef ImageContainerBySTLVector<Domain, float> FloatImage;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef DistanceToMeasureEdge<FloatImage> Distance;

    typedef ImageContainerByITKImage<Domain, int> ITKImage;

    typedef DGtal::HessianRecursiveGaussian<GrayLevelImage> MyHessian;
    typedef typename MyHessian::OutputImage HessianImage;
    typedef ImageToDistance<HessianImage> MyDistanceImage;
    typedef DGtal::FrangiVesselness<HessianImage> Vesselness;

    GrayLevelImage img = DGtal::ITKReader<GrayLevelImage>::importITK(inputFilename);
    auto domain = img.domain();
    FloatImage fimg(img.domain());
    FloatImage::Iterator outIt = fimg.begin();
    for (GrayLevelImage::ConstIterator it = img.begin(), itE = img.end(); it != itE; ++it) {
        float v = ((float) *it) * 1.0f / thresholdMax;
        *outIt++ = v;
    }
    trace.beginBlock("Computing delta-distance.");
    Distance delta;
    if (vm.count("distance")) {
        string distanceName = vm["distance"].as<std::string>();
        FloatImage distanceImage = ITKReader<FloatImage>::importITK(distanceName);
        FloatImage distanceImage2(distanceImage.domain());
        for (const Point& p : distanceImage.domain()) {
            float value = distanceImage(p);
            distanceImage2.setValue(p, value*value);
        }
        delta = Distance(distanceImage2);

    }
    else {
        delta = Distance(mass, fimg, rmax, 0.0);
    }

    trace.endBlock();

    std::map<Point, double> pToValues;
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();


    ITKImage image = DGtal::ITKReader<ITKImage>::importITK(inputFilename);



    double maxDistance = delta(
            *std::max_element(domain.begin(), domain.end(), [&](const Point &p1, const Point &p2) {
                return delta.distance(p1) < delta.distance(p2);
            }));

    std::vector<MyDistanceImage> hessianFiltersVector;

    DGtal::trace.info() << "Max distance= " << maxDistance << std::endl;

    DGtal::trace.beginBlock("Multiscale hessian");
    MyHessian hessian(img, 1.0, false);
    HessianImage hessianImage(img.domain());
    for (float i = 0.5; i <= maxDistance; i+=0.5) {
        DGtal::trace.progressBar(i, maxDistance);
        double currentSigma = i / 2.f;
        hessian.setSigma(currentSigma);
        auto res = hessian.computeHessian();
        for (const auto& p : image.domain()) {
                const double distanceP = delta.distance(p);
                if (distanceP >= i && distanceP < i + 0.5)
                        hessianImage.setValue(p, res(p));
        }
    }
    DGtal::trace.endBlock();


    DGtal::trace.beginBlock("Frangi vesselness");
    Vesselness vesselness(hessianImage);
    auto vesselnessImage = vesselness.computeVesselness();
    DGtal::trace.endBlock();

//    GradientColorMap<float> cmap_grad(min, max);
//
//    cmap_grad.addColor(Color(255, 255, 255));
//    cmap_grad.addColor(Color(255, 255, 0));
//    cmap_grad.addColor(Color(255, 0, 0));
//    cmap_grad.addColor(Color(0, 255, 0));
//    cmap_grad.addColor(Color(0, 0, 255));
//    cmap_grad.addColor(Color(0, 0, 0));


//    for (const auto &pair : pToValues) {
//        Point p = pair.first;
//        double value = pair.second;
//        Z3i::Point p3D;
//        for (int i = 0; i < dimension; i++) {
//            p3D[i] = p[i];
//        }
//        Color currentColor = cmap_grad(value);
//        currentColor.alpha(value * 255);
//        viewer << CustomColors3D(currentColor, currentColor) << p3D;
//    }



    size_t lastindex = outname.find_last_of(".");
    string extension = outname.substr(outname.find_last_of("."));
    string outRadiusName = outname.substr(0, lastindex) + "_radius" + extension;

    FloatImage out(vesselnessImage.domain());
    for (const Z3i::Point &p : out.domain()) {
        out.setValue(p, (float) vesselnessImage(p));
    }
    DGtal::ITKWriter<FloatImage>::exportITK(outname, out);
    app.exec();

    return 0;
}
