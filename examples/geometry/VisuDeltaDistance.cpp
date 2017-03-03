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

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "geometry/DistanceToMeasureEdge.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main(int argc, char **argv) {
    using namespace DGtal;
    using namespace DGtal::Z2i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage2D;
    typedef ImageContainerBySTLVector<Domain, float> FloatImage2D;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef DistanceToMeasureEdge<FloatImage2D> Distance;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
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
    int thresholdMax = vm["thresholdMax"].as<int>();
    double mass = vm["mass"].as<double>();
    double rmax = vm["rmax"].as<double>();


    GrayLevelImage2D img = GenericReader<GrayLevelImage2D>::import(inputFilename);
    auto domain = img.domain();
    FloatImage2D fimg(img.domain());
    FloatImage2D::Iterator outIt = fimg.begin();
    for (GrayLevelImage2D::ConstIterator it = img.begin(), itE = img.end(); it != itE; ++it) {
        float v = ((float) *it) * 1.0f / thresholdMax;
        *outIt++ = v;
    }
    trace.beginBlock("Computing delta-distance.");
    Distance delta(mass, fimg, rmax);
    trace.endBlock();


    QApplication app(argc, argv);
    ViewerDistanceBall<Distance> viewer(delta);
    viewer.show();

    double min = fimg(*std::min_element(domain.begin(), domain.end(), [&](const Z2i::Point &p1,
                                                                          const Z2i::Point &p2) {
        return fimg(p1) < fimg(p2);
    }));

    double max = fimg(*std::max_element(domain.begin(), domain.end(), [&](const Z2i::Point &p1,
                                                                          const Z2i::Point &p2) {
        return fimg(p1) < fimg(p2);
    }));
    DGtal::trace.info() << min << " " << max << std::endl;
    GradientColorMap<float> cmap_grad(min, max);

    cmap_grad.addColor(Color(255, 255, 255));
    cmap_grad.addColor(Color(255, 255, 0));
    cmap_grad.addColor(Color(255, 0, 0));
    cmap_grad.addColor(Color(0, 255, 0));
    cmap_grad.addColor(Color(0, 0, 255));
    cmap_grad.addColor(Color(0, 0, 0));


    for (const Point &p : delta.domain()) {
        Z3i::Point p3D(p[0], p[1], 0);
        double value = fimg(p);
        Color currentColor = cmap_grad(value);
        viewer << CustomColors3D(currentColor, currentColor) << p3D;
    }


    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
