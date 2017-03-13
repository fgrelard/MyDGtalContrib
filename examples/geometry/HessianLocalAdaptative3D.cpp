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


int main(int argc, char **argv) {
    using namespace DGtal;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
            ("output,o", po::value<std::string>(), "vol file (corresponding volume)")
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
    typedef int InputPixelType;
    typedef double OutputPixelType;
    typedef itk::Image<InputPixelType, dimension> InputImageType;
    typedef itk::Image<OutputPixelType, dimension> OutputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageRegionIteratorWithIndex<OutputImageType> ImageIterator;
    typedef itk::DiscreteHessianGaussianImageFunction<InputImageType>
            HessianFilterType;
    typedef itk::Point<InputPixelType, dimension> ITKPoint;
    typedef typename InputImageType::IndexType IndexType;
    typedef ImageContainerByITKImage<Domain, int> ITKImage;
    typedef HessianFilterType::TensorType::EigenValuesArrayType EigenValues;
    typedef itk::SymmetricSecondRankTensor<double, dimension> HessianPixelType;
    typedef itk::Image<HessianPixelType, dimension> HessianImageType;


    GrayLevelImage img = DGtal::ITKReader<GrayLevelImage>::importITK(inputFilename);
    GrayLevelImage imgInvert = invert(img);
    auto domain = img.domain();
    FloatImage fimg(imgInvert.domain());
    FloatImage::Iterator outIt = fimg.begin();
    for (GrayLevelImage::ConstIterator it = imgInvert.begin(), itE = imgInvert.end(); it != itE; ++it) {
        float v = ((float) *it) * 1.0f / thresholdMax;
        *outIt++ = v;
    }
    trace.beginBlock("Computing delta-distance.");
    Distance delta(mass, fimg, rmax);
    trace.endBlock();

    std::map<Point, double> pToValues;
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();


    ITKImage image = DGtal::ITKReader<ITKImage>::importITK(inputFilename);
    ITKImage::ITKImagePointer imagePointer = image.getITKImagePointer();
    HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
    HessianImageType::Pointer hessianImage = HessianImageType::New();
    HessianImageType::RegionType region;
    HessianImageType::IndexType start;
    for (int i = 0; i < dimension; i++)
        start[i] = 0;

    HessianImageType::SizeType size;
    for (int i = 0; i < dimension; i++) {
        size[i] = imagePointer->GetRequestedRegion().GetSize()[i];
    }
    region.SetSize(size);
    region.SetIndex(start);
    hessianImage->SetRegions(region);
    hessianImage->Allocate();
    hessianImage->Update();

    double maxDistance = delta.distance(
            *std::max_element(domain.begin(), domain.end(), [&](const Point &p1, const Point &p2) {
                return delta.distance(p1) < delta.distance(p2);
            }));

    std::vector<HessianFilterType::Pointer> hessianFiltersVector;
    for (int i = 0; i <= maxDistance; i++) {
        HessianFilterType::Pointer hessianFilter = HessianFilterType::New();

        hessianFilter->SetUseImageSpacing(false);
        hessianFilter->SetInputImage(imagePointer);
        hessianFilter->SetNormalizeAcrossScale(true);
        hessianFilter->SetSigma(i + 1);
        hessianFilter->Initialize();
        hessianFiltersVector.push_back(hessianFilter);
    }


    for (const auto &p : image.domain()) {
        EigenValues eigenValues;
        ITKPoint itkP;
        HessianImageType::IndexType index;

        for (int i = 0; i < dimension; i++) {
            itkP[i] = p[i];
            index[i] = p[i];
        }
        const double distanceP = delta.distance(p);
        HessianFilterType::Pointer hessianFilter = hessianFiltersVector[(int) distanceP];
        auto hessian = hessianFilter->Evaluate(itkP);
        hessian.ComputeEigenValues(eigenValues);
        hessianImage->SetPixel(index, hessian);
    }

    typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, OutputImageType>
            ObjectnessFilterType;
    ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
    objectnessFilter->SetBrightObject(true);
    objectnessFilter->SetScaleObjectnessMeasure(false);
    objectnessFilter->SetAlpha(0.5);
    objectnessFilter->SetBeta(1.0);
    objectnessFilter->SetGamma(5.0);
    objectnessFilter->SetInput(hessianImage);
    objectnessFilter->Update();
    OutputImageType::Pointer outObject = objectnessFilter->GetOutput();

    ImageIterator it(outObject, outObject->GetRequestedRegion());
    it.GoToBegin();
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    while (!it.IsAtEnd()) {
        double value = it.Get();
        if (value < min)
            min = value;
        if (value > max)
            max = value;
        OutputImageType::IndexType itkP = it.GetIndex();
        Point p;
        for (int i = 0; i < dimension; i++) {
            p[i] = itkP[i];
        }
        pToValues[p] = value;
        ++it;
    }
    DGtal::trace.info() << min << " " << max << std::endl;
    GradientColorMap<float> cmap_grad(min, max);

    cmap_grad.addColor(Color(255, 255, 255));
    cmap_grad.addColor(Color(255, 255, 0));
    cmap_grad.addColor(Color(255, 0, 0));
    cmap_grad.addColor(Color(0, 255, 0));
    cmap_grad.addColor(Color(0, 0, 255));
    cmap_grad.addColor(Color(0, 0, 0));


    for (const auto &pair : pToValues) {
        Point p = pair.first;
        double value = pair.second;
        Z3i::Point p3D;
        for (int i = 0; i < dimension; i++) {
            p3D[i] = p[i];
        }
        Color currentColor = cmap_grad(value);
        currentColor.alpha(value * 255);
        viewer << CustomColors3D(currentColor, currentColor) << p3D;
    }

    GrayLevelImage out(img.domain());
    for (const auto &pair : pToValues) {
        Point p = pair.first;
        unsigned char value = pair.second * 255 / max;
        out.setValue(p, value);
    }


    double maxRadius = delta.distance(*std::max_element(domain.begin(), domain.end(), [&](const Point &one,
                                                                                          const Point &two) {
        return delta.distance(one) < delta.distance(two);
    }));

    GrayLevelImage outRadius(img.domain());
    for (const Point &p : domain) {
        unsigned char value = delta.distance(p) * 255 / maxRadius;
        outRadius.setValue(p, value);
    }

    ITKWriter<GrayLevelImage>::exportITK(outname, out);
    ITKWriter<GrayLevelImage>::exportITK(outRadiusName, outRadius);
    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}

