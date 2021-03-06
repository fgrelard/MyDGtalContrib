#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/geometry/volumes/distance/VoronoiMap.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkDiscreteHessianGaussianImageFunction.h>
#include <DGtal/io/writers/ITKWriter.h>
#include <itkHessianToObjectnessMeasureImageFilter.h>
#include <itkMultiScaleHessianBasedMeasureImageFilter.h>
#include <itkHessian3DToVesselnessMeasureImageFilter.h>
#include <itktools/itkModifiedKrissianVesselnessImageFilter.h>
#include <itktools/itkMultiScaleGaussianEnhancementImageFilter.h>
#include <itkImageDuplicator.h>
#include <QtWidgets/QApplication>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/viewers/Viewer3D.h"
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

struct DoubleToFloatFunctor {

    float operator()(double value) {
        return (float) value;
    }
};


int main(int argc, char **argv) {
    using namespace DGtal;
    using namespace DGtal::Z3i;

    typedef ImageContainerBySTLVector<Domain, double> GrayLevelImage2D;

    typedef ImageContainerBySTLVector<Domain, double> FloatImage2D;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef double InputPixelType;
    typedef double OutputPixelType;
    typedef itk::Image<InputPixelType, 3> InputImageType;
    typedef itk::Image<OutputPixelType, 3> OutputImageType;

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageRegionIteratorWithIndex<OutputImageType> ImageIterator;
    typedef itk::Point<InputPixelType, 3> ITKPoint;
    typedef typename InputImageType::IndexType IndexType;
    typedef ImageContainerByITKImage<Domain, double> ITKImage;
    typedef ImageContainerByITKImage<Domain, double> ITKImageDouble;
    typedef ImageContainerByITKImage<Domain, float> ITKImageFloat;
    typedef itk::SymmetricSecondRankTensor<double, 3> HessianPixelType;
    typedef itk::Image<HessianPixelType, 3> HessianImageType;
    typedef itk::MultiScaleHessianBasedMeasureImageFilter<InputImageType, HessianImageType, OutputImageType>
            MultiVesselnessType;

    //Frangi
    typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, OutputImageType> ObjectnessFilterType;

    // Sato
//    typedef itk::Hessian3DToVesselnessMeasureImageFilter< OutputPixelType >
//            ObjectnessFilterType;

    // Krissian
//    typedef itk::ModifiedKrissianVesselnessImageFilter<HessianImageType, OutputImageType>
//            ObjectnessFilterType;
//    typedef itk::MultiScaleGaussianEnhancementImageFilter<InputImageType, OutputImageType>
//            MultiVesselnessType;

    typedef typename MultiVesselnessType::ScalesImageType ScaleImageType;
//    typedef typename MultiVesselnessType::EigenValueArrayType EigenValueArrayType;
//    typedef itk::Functor::ModifiedKrissianVesselnessFunctor<EigenValueArrayType, OutputPixelType> KrissianFunctor;

    typedef itk::HessianRecursiveGaussianImageFilter<InputImageType, HessianImageType> HessianFilter;


    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
            ("output,o", po::value<std::string>(), "vol file (corresponding volume)")
        ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
         ("numberSteps,n", po::value<int>()->default_value(10), "maximum threshold for binarization")
        ("sigmaMin,s", po::value<double>()->default_value(1.0), "minimum sigma for Frangi")
        ("sigmaMax,S", po::value<double>()->default_value(3.0), "maximum sigma for Frangi")
        ("alpha,a", po::value<double>()->default_value(0.1), "Alpha")
        ("beta,b", po::value<double>()->default_value(1.0), "Beta")
        ("gamma,g", po::value<double>()->default_value(10.0), "Gamma");


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
    double alpha = vm["alpha"].as<double>();
    double beta = vm["beta"].as<double>();
    double gamma = vm["gamma"].as<double>();
    double sigmaMin = vm["sigmaMin"].as<double>();
    double sigmaMax = vm["sigmaMax"].as<double>();
    double numberSteps = vm["numberSteps"].as<int>();

    GrayLevelImage2D img = DGtal::ITKReader<GrayLevelImage2D>::importITK(inputFilename);
    auto domain = img.domain();
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();

    ITKImage image = DGtal::ITKReader<ITKImage>::importITK(inputFilename);
    ITKImage::ITKImagePointer imagePointer = image.getITKImagePointer();
    size_t lastindex = outname.find_last_of(".");
    string extension = outname.substr(outname.find_last_of("."));

    ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();

    // Frangi
    objectnessFilter->SetBrightObject(true);
    objectnessFilter->SetScaleObjectnessMeasure(false);
    objectnessFilter->SetAlpha(alpha);
    objectnessFilter->SetBeta(beta);
    objectnessFilter->SetGamma(gamma);

    // Sato
//    objectnessFilter->SetAlpha1(0.5);
//    objectnessFilter->SetAlpha2(2.0);

    //Krissian
//    typename KrissianFunctor::Pointer functor = KrissianFunctor::New();
//    functor->SetBrightObject(true);

    MultiVesselnessType::Pointer vesselnessType = MultiVesselnessType::New();

    // Krissian
//    vesselnessType->SetUnaryFunctor(functor);

    //Sato+Frangi
    vesselnessType->SetHessianToMeasureFilter(objectnessFilter);

    //All
    vesselnessType->SetInput(imagePointer);
    vesselnessType->SetNonNegativeHessianBasedMeasure(false);
    vesselnessType->SetGenerateScalesOutput(true);
    vesselnessType->SetSigmaMinimum(sigmaMin);
    vesselnessType->SetSigmaMaximum(sigmaMax);
    vesselnessType->SetNumberOfSigmaSteps(numberSteps);
    vesselnessType->PrepareOutputs();
    vesselnessType->UpdateLargestPossibleRegion();
    vesselnessType->UpdateOutputInformation();
    vesselnessType->Update();
    OutputImageType::Pointer outObject = vesselnessType->GetOutput();
    ITKImageDouble outDouble(outObject);
    ITKImageFloat outFloat(outDouble.domain());
    for (const Point &p : outDouble.domain()) {
        outFloat.setValue(p, (float) outDouble(p));
    }
    string outRadiusName = outname.substr(0, lastindex) + "_unique" + extension;
    ITKWriter<ITKImageFloat>::exportITK(outRadiusName, outFloat);


    ScaleImageType::ConstPointer outScaleObject = vesselnessType->GetScalesOutput();
    typedef itk::ImageDuplicator<ScaleImageType> DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(outScaleObject);
    duplicator->Update();
    ScaleImageType::Pointer copy = duplicator->GetOutput();

    //Frangi+Sato
    ITKImageFloat outScaleFloat(copy);

    // Krissian
//    ITKImageDouble outScaleDouble(copy);
//    ITKImageFloat outScaleFloat(outScaleDouble.domain());
//    for (const Point &p : outScaleDouble.domain()) {
//        outScaleFloat.setValue(p, (float) outScaleDouble(p));
//    }
    string outScaleName = outname.substr(0, lastindex) + "_scale" + extension;
    ITKWriter<ITKImageFloat>::exportITK(outScaleName, outScaleFloat);


    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
