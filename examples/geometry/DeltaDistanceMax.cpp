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
#include <geometry/DistanceToMeasureRelevantScale.h>
#include <viewer/ViewerDistanceBall.h>
#include <DGtal/io/writers/ITKWriter.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "geometry/DistanceToMeasureEdge.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

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

    ImageToDistance() : myDistance(std::numeric_limits<float>::max()) {}

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
    double valueMax = std::numeric_limits<double>::max();
    for (const Value& imageToDistance : vImage) {
        float diff = std::abs(imageToDistance.myDistance - distance);
        if (diff < valueMax) {
            valueMax = diff;
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
        ("distance,d", po::value<std::string>(), "distance image")
        ("output,o", po::value<std::string>(), "image file (corresponding volume)")
        ("mass,a", po::value<double>()->default_value(1), "Mass to integrate for distance to measure")
        ("alpha,f", po::value<double>()->default_value(0.1), "Alpha")
        ("beta,s", po::value<double>()->default_value(1.0), "Beta")
        ("gamma,t", po::value<double>()->default_value(10.0), "Gamma")
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


    string outname = vm["output"].as<std::string>();

    int thresholdMax = vm["thresholdMax"].as<int>();
    double mass = vm["mass"].as<double>();
    double rmax = vm["rmax"].as<double>();
    double alpha = vm["alpha"].as<double>();
    double beta = vm["beta"].as<double>();
    double gamma = vm["gamma"].as<double>();
    constexpr unsigned int dimension = 3;

    typedef SpaceND<dimension> Space;
    typedef HyperRectDomain<Space> Domain;
    typedef typename Space::Point Point;
    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage;

    typedef ImageContainerBySTLVector<Domain, float> FloatImage;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef DistanceToMeasureRelevantScale<FloatImage> Distance;
    typedef int InputPixelType;
    typedef double OutputPixelType;

    typedef ImageContainerByITKImage<Domain, int> ITKImage;


    trace.beginBlock("Computing delta-distance.");
    string distanceName = vm["distance"].as<std::string>();
    FloatImage distanceImage = ITKReader<FloatImage>::importITK(distanceName);
    FloatImage distanceImage2(distanceImage.domain());
    for (const Point& p : distanceImage.domain()) {
        float value = distanceImage(p);
        distanceImage2.setValue(p, value*value);
    }
    Distance delta(distanceImage2);
    trace.endBlock();

    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();

    FloatImage distanceImageMax(distanceImage.domain());
    int i = 0;
     for (const Point& p : distanceImage.domain()) {
         trace.progressBar(i, distanceImage.domain().size());
         float valueMax = delta.maxAlongProj(p);
         distanceImageMax.setValue(p, valueMax);
         i++;
    }


    ITKWriter<FloatImage>::exportITK(outname, distanceImageMax);
    //app.exec();

    return 0;
}
