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
    typedef unsigned char InputPixelType;
    typedef double OutputPixelType;
    typedef itk::Image<InputPixelType, 2> InputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageRegionIteratorWithIndex<InputImageType> ImageIterator;
    typedef itk::DiscreteHessianGaussianImageFunction<InputImageType>
            HessianFilterType;
    typedef itk::Point<InputPixelType, 2> ITKPoint;
    typedef typename InputImageType::IndexType IndexType;
    typedef ImageContainerByITKImage<Domain, unsigned char> ITKImage;

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


    GrayLevelImage2D img = DGtal::ITKReader<GrayLevelImage2D>::importITK(inputFilename);
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

    std::map<Point, double> pToValues;
    QApplication app(argc, argv);
    ViewerDistanceBall<Distance> viewer(delta);
    viewer.show();


    ITKImage image = DGtal::ITKReader<ITKImage>::importITK(inputFilename);
    ITKImage::ITKImagePointer imagePointer = image.getITKImagePointer();

    HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
    hessianFilter->SetInputImage(imagePointer);
    ImageIterator it(imagePointer, imagePointer->GetRequestedRegion());
    hessianFilter->SetSigma(1.0);
    hessianFilter->Initialize();

    for (const auto &p : image.domain()) {
        ITKPoint itkP;
        itkP = p[0];
        itkP = p[1];
//        IndexType index = it.GetIndex();
//        imagePointer->TransformIndexToPhysicalPoint(index, itkP);
//        DGtal::trace.info() << index << std::endl;
        //Z2i::Point p(itkP[0], itkP[1]);
        const double distanceP = delta.distance(p);

        auto hessian = hessianFilter->Evaluate(itkP);
        double max = *std::max_element(hessian.Begin(), hessian.End(), [&](const double &first, const double &second) {
            return first * first < second * second;
        });
        double max2 = max * max;
        pToValues[p] = max2;
    }

    double min = std::min_element(pToValues.begin(), pToValues.end(), [&](const std::pair<Z2i::Point, double> &p1,
                                                                          const std::pair<Z2i::Point, double> &p2) {
        return p1.second < p2.second;
    })->second;

    double max = std::max_element(pToValues.begin(), pToValues.end(), [&](const std::pair<Z2i::Point, double> &p1,
                                                                          const std::pair<Z2i::Point, double> &p2) {
        return p1.second < p2.second;
    })->second;
    DGtal::trace.info() << min << " " << max << std::endl;
    GradientColorMap<float> cmap_grad(min, max);

    cmap_grad.addColor(Color(255, 255, 255));
    cmap_grad.addColor(Color(255, 255, 0));
    cmap_grad.addColor(Color(255, 0, 0));
    cmap_grad.addColor(Color(0, 255, 0));
    cmap_grad.addColor(Color(0, 0, 255));
    cmap_grad.addColor(Color(0, 0, 0));


    for (const auto &pair : pToValues) {
        Z2i::Point p = pair.first;
        double value = pair.second;
        Z3i::Point p3D(p[0], p[1], 0);
        Color currentColor = cmap_grad(value);
        viewer << CustomColors3D(currentColor, currentColor) << p3D;
    }


    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
