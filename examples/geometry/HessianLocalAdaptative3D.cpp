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
        ("input,i", po::value<std::string>(), "image file (corresponding volume)")
        ("distance,d", po::value<std::string>(), "distance image")
        ("localDistance,l", po::value<std::string>(), "distance image")
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
    typedef DistanceToMeasureEdge<FloatImage> Distance;
    typedef int InputPixelType;
    typedef double OutputPixelType;
    typedef itk::Image<InputPixelType, dimension> InputImageType;
    typedef itk::Image<OutputPixelType, dimension> OutputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageRegionIteratorWithIndex<OutputImageType> ImageIterator;
    typedef itk::SymmetricSecondRankTensor<double, dimension> HessianPixelType;
    typedef itk::Image<HessianPixelType, dimension> HessianImageType;

    typedef itk::HessianRecursiveGaussianImageFilter<InputImageType> HessianFilterType;

    typedef itk::Point<InputPixelType, dimension> ITKPoint;
    typedef typename InputImageType::IndexType IndexType;
    typedef ImageContainerByITKImage<Domain, int> ITKImage;
    typedef itk::ImageRegionIteratorWithIndex<HessianImageType> HessianIterator;
    typedef typename HessianImageType::IndexType HessianIndexType;
    typedef ImageToDistance<HessianImageType::Pointer> DistanceImage;

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
    double maxDistance = delta(
        *std::max_element(domain.begin(), domain.end(), [&](const Point &p1, const Point &p2) {
                return delta.distance(p1) < delta.distance(p2);
            }));
    FloatImage localDistanceImage(domain);
    if (vm.count("localDistance")) {
        std::string localDistanceName = vm["localDistance"].as<std::string>();
        FloatImage localDistanceImage = ITKReader<FloatImage>::importITK(localDistanceName);
    }
    else {
        for (const Z3i::Point&  p : domain)
            localDistanceImage.setValue(p, maxDistance);
    }

    trace.endBlock();

    std::map<Point, double> pToValues;
    //    QApplication app(argc, argv);
    //    Viewer3D<> viewer;
    //    viewer.show();


    ITKImage image = DGtal::ITKReader<ITKImage>::importITK(inputFilename);

    ITKImage::ITKImagePointer imagePointer = image.getITKImagePointer();
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


//    double maxDistance = 1.0;
    std::vector<DistanceImage> hessianFiltersVector;

    DGtal::trace.info() << "Max distance= " << maxDistance << std::endl;

    DGtal::trace.beginBlock("Multiscale hessian");
    HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
    hessianFilter->SetInput(imagePointer);
    hessianFilter->SetNormalizeAcrossScale(true);
    for (float i = 1; i <= maxDistance-1; i++) {
        DGtal::trace.progressBar(i, maxDistance);
//        double currentSigma = (maxDistance * i)/(2 * (i + 0.5));
        double currentSigma = i/(maxDistance-i);
//        currentSigma = i /2.0;
        DGtal::trace.info() << currentSigma << " " << i << std::endl;
        hessianFilter->SetSigma(currentSigma);
        hessianFilter->Update();
        HessianImageType::Pointer currentHessianImage = hessianFilter->GetOutput();
        DGtal::trace.beginBlock("Hessian computation");
        int j = 0;

        for (const auto &p : image.domain()) {
            if (delta.distance(p) == 0.0) continue;
            trace.progressBar(i, image.domain().size());
            j++;

            ITKPoint itkP;
            HessianImageType::IndexType index;

            for (int i = 0; i < dimension; i++) {
                index[i] = p[i];
            }
            double distanceP;
            if (localDistanceImage(p) > 0)
                 distanceP = delta.distance(p) * (delta.distance(p) / localDistanceImage(p));
            else
                distanceP = delta.distance(p);
            size_t value = std::round(distanceP);
            if (value == i) {
                HessianPixelType hessian = currentHessianImage->GetPixel(index);
                hessianImage->SetPixel(index, hessian);
            }
        }
        DGtal::trace.endBlock();

        // DistanceImage imageToDistance(hessianImage, i);
        // hessianFiltersVector.push_back(imageToDistance);
        // hessianImage->DisconnectPipeline();
    }
    DGtal::trace.endBlock();



    DGtal::trace.beginBlock("Frangi vesselness");
    typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, OutputImageType>
            ObjectnessFilterType;
    ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
    objectnessFilter->SetBrightObject(true);
    objectnessFilter->SetScaleObjectnessMeasure(false);
    objectnessFilter->SetAlpha(alpha);
    objectnessFilter->SetBeta(beta);
    objectnessFilter->SetGamma(gamma);
    objectnessFilter->SetObjectDimension(1);
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
    DGtal::trace.endBlock();

    DGtal::trace.info() << min << " " << max << std::endl;
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

    FloatImage out(img.domain());
    for (const auto &pair : pToValues) {
        Point p = pair.first;
        unsigned char value = pair.second * 255 / max;
        out.setValue(p, pair.second);
    }


    size_t lastindex = outname.find_last_of(".");
    string extension = outname.substr(outname.find_last_of("."));
    string outRadiusName = outname.substr(0, lastindex) + "_radius" + extension;
    FloatImage outRadius(img.domain());
    for (const Point &p : domain) {
        float value = delta.distance(p);
        //float value = 1.0;
        outRadius.setValue(p, value);
    }

    ITKWriter<FloatImage>::exportITK(outname, out);
    ITKWriter<FloatImage>::exportITK(outRadiusName, outRadius);
    //app.exec();

    return 0;
}

