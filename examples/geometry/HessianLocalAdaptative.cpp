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
#include <DGtal/io/writers/ITKWriter.h>
#include <itkHessianToObjectnessMeasureImageFilter.h>
#include <hessian/DiscreteHessianFunction.h>
#include <itkDiscreteHessianGaussianImageFunction.h>
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
    using namespace DGtal::Z2i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage2D;

    typedef ImageContainerBySTLVector<Domain, float> FloatImage2D;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef DistanceToMeasureEdge<FloatImage2D> Distance;
    typedef int InputPixelType;
    typedef double OutputPixelType;
    typedef itk::Image<InputPixelType, 2> InputImageType;
    typedef itk::Image<OutputPixelType, 2> OutputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageRegionIteratorWithIndex<OutputImageType> ImageIterator;
    typedef itk::DiscreteHessianGaussianImageFunction<InputImageType>
            HessianFilterType;
    typedef itk::Point<InputPixelType, 2> ITKPoint;
    typedef typename InputImageType::IndexType IndexType;
    typedef ImageContainerByITKImage<Domain, int> ITKImage;
    typedef ImageContainerByITKImage<Domain, double> ITKImageDouble;

    typedef HessianFilterType::TensorType::EigenValuesArrayType EigenValues;
    typedef itk::SymmetricSecondRankTensor<double, 2> HessianPixelType;
    typedef itk::Image<HessianPixelType, 2> HessianImageType;
    typedef DGtal::DiscreteHessianFunction<ITKImage> DiscreteHessian;

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


    GrayLevelImage2D img = DGtal::ITKReader<GrayLevelImage2D>::importITK(inputFilename);
    GrayLevelImage2D imgInvert = invert(img);
    auto domain = img.domain();
    FloatImage2D fimg(imgInvert.domain());
    FloatImage2D::Iterator outIt = fimg.begin();
    for (GrayLevelImage2D::ConstIterator it = imgInvert.begin(), itE = imgInvert.end(); it != itE; ++it) {
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
    ITKPoint origin;
    origin[0] = 0;
    origin[1] = 0;
    ITKImage::ITKImagePointer imagePointer = image.getITKImagePointer();
    HessianImageType::Pointer hessianImage = HessianImageType::New();
    HessianImageType::RegionType region;
    HessianImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;

    HessianImageType::SizeType size;
    size[0] = imagePointer->GetRequestedRegion().GetSize()[0];
    size[1] = imagePointer->GetRequestedRegion().GetSize()[1];

    region.SetSize(size);
    region.SetIndex(start);
    hessianImage->SetRegions(region);
    hessianImage->Allocate();
    hessianImage->Update();

    double maxDistance = delta.distance(
            *std::max_element(domain.begin(), domain.end(), [&](const Point &p1, const Point &p2) {
                return delta.distance(p1) < delta.distance(p2);
            }));

    maxDistance = 0;
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

    DiscreteHessian hessian(image, 1);
    typename DiscreteHessian::OutputImage matrix = hessian.computeHessian();

    for (const auto &p : image.domain()) {
        EigenValues eigenValues;
        ITKPoint itkP;
        itkP[0] = p[0];
        itkP[1] = p[1];
        double distanceP = delta.distance(p);
        distanceP = 0;
        HessianFilterType::Pointer hessianFilter = hessianFiltersVector[(int) distanceP];
//        hessianFilter->SetSigma(distanceP);
//        hessianFilter->Initialize();
        auto hessian = hessianFilter->Evaluate(itkP);
        hessian.ComputeEigenValues(eigenValues);
        HessianImageType::IndexType index;
        index[0] = p[0];
        index[1] = p[1];
        auto value = matrix(p);
        DGtal::trace.info() << "ITK= ";
        for (int i = 0; i < 3; i++) {
            DGtal::trace.info() << hessian[i] << " ";
        }
        DGtal::trace.info() << std::endl;

//        DGtal::trace.info() << "DGtal= ";
//        for (const double &d : value) {
//            DGtal::trace.info() << d << " ";
//        }
        DGtal::trace.info() << std::endl;
//        double *a = &value[0];
//        hessianImage->SetPixel(index, a);
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
        Point p(itkP[0], itkP[1]);
        pToValues[p] = value;
        ++it;
    }
//    double min = std::min_element(pToValues.begin(), pToValues.end(), [&](const std::pair<Z2i::Point, double> &p1,
//                                                                          const std::pair<Z2i::Point, double> &p2) {
//        return p1.second < p2.second;
//    })->second;
//
//    double max = std::max_element(pToValues.begin(), pToValues.end(), [&](const std::pair<Z2i::Point, double> &p1,
//                                                                          const std::pair<Z2i::Point, double> &p2) {
//        return p1.second < p2.second;
//    })->second;
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

    GrayLevelImage2D out(img.domain());
    for (const auto &pair : pToValues) {
        Z2i::Point p = pair.first;
        unsigned char value = pair.second * 255 / max;
        out.setValue(p, value);
    }


    ITKWriter<GrayLevelImage2D>::exportITK(outname, out);

    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
