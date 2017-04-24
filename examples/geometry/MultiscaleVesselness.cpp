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
#include <itkHessian3DToVesselnessMeasureImageFilter.h>
#include <itktools/itkModifiedKrissianVesselnessImageFilter.h>
#include <itktools/itkMultiScaleGaussianEnhancementImageFilter.h>
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

struct DoubleToFloatFunctor {

    float operator()(double value) {
        return (float) value;
    }
};


int main(int argc, char **argv) {
    using namespace DGtal;
    using namespace DGtal::Z3i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage2D;

    typedef ImageContainerBySTLVector<Domain, double> FloatImage2D;
    typedef ImageContainerBySTLVector<Domain, DGtal::Color> OutImage;
    typedef DistanceToMeasureEdge<FloatImage2D> Distance;
    typedef int InputPixelType;
    typedef double OutputPixelType;
    typedef itk::Image<InputPixelType, 3> InputImageType;
    typedef itk::Image<OutputPixelType, 3> OutputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageRegionIteratorWithIndex<OutputImageType> ImageIterator;
    typedef itk::Point<InputPixelType, 3> ITKPoint;
    typedef typename InputImageType::IndexType IndexType;
    typedef ImageContainerByITKImage<Domain, int> ITKImage;
    typedef ImageContainerByITKImage<Domain, double> ITKImageDouble;
    typedef ImageContainerByITKImage<Domain, float> ITKImageFloat;
    typedef itk::SymmetricSecondRankTensor<double, 3> HessianPixelType;
    typedef itk::Image<HessianPixelType, 3> HessianImageType;

    //Frangi
    typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, OutputImageType> ObjectnessFilterType;

    // Sato
//    typedef itk::Hessian3DToVesselnessMeasureImageFilter<OutputPixelType>
//            ObjectnessFilterType;

    typedef itk::MultiScaleHessianBasedMeasureImageFilter<InputImageType, HessianImageType, OutputImageType>
            MultiVesselnessType;

    // Krissian
//    typedef itk::ModifiedKrissianVesselnessImageFilter< HessianImageType, OutputImageType >
//            ObjectnessFilterType;
//    typedef itk::MultiScaleGaussianEnhancementImageFilter<InputImageType, OutputImageType>
//            MultiVesselnessType;
//    typedef typename MultiVesselnessType::EigenValueArrayType                             EigenValueArrayType;
//    typedef itk::Functor::ModifiedKrissianVesselnessFunctor<EigenValueArrayType, OutputPixelType > KrissianFunctor;

    typedef itk::HessianRecursiveGaussianImageFilter<InputImageType, HessianImageType> HessianFilter;


    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
            ("output,o", po::value<std::string>(), "vol file (corresponding volume)")
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


    GrayLevelImage2D img = DGtal::ITKReader<GrayLevelImage2D>::importITK(inputFilename);
    auto domain = img.domain();
    std::map<Point, double> pToValues;
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();

    ITKImage image = DGtal::ITKReader<ITKImage>::importITK(inputFilename);
    ITKImage::ITKImagePointer imagePointer = image.getITKImagePointer();
    size_t lastindex = outname.find_last_of(".");
    string extension = outname.substr(outname.find_last_of("."));

    float maxSigma = 5.0;
    ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();

    DGtal::trace.beginBlock("Multiscale hessian");
    for (float i = 0.1; i <= maxSigma; i += 0.1) {
        DGtal::trace.progressBar(i, maxSigma);

        // Frangi
        objectnessFilter->SetBrightObject(true);
        objectnessFilter->SetScaleObjectnessMeasure(false);
        objectnessFilter->SetAlpha(0.5);
        objectnessFilter->SetBeta(1.0);
        objectnessFilter->SetGamma(5.0);

        // Sato
//        objectnessFilter->SetAlpha1(0.5);
//        objectnessFilter->SetAlpha2(2.0);


        //Krissian
//        typename KrissianFunctor::Pointer functor = KrissianFunctor::New();
//        functor->SetBrightObject( true );

        MultiVesselnessType::Pointer vesselnessType = MultiVesselnessType::New();

        //Krissian
        //vesselnessType->SetUnaryFunctor( functor );

        //Frangi+Sato
        vesselnessType->SetHessianToMeasureFilter(objectnessFilter);


        vesselnessType->SetInput(imagePointer);
        vesselnessType->SetNonNegativeHessianBasedMeasure(true);
        vesselnessType->SetSigmaMaximum(i);
        vesselnessType->SetSigmaMinimum(i);
        vesselnessType->SetNumberOfSigmaSteps(1);
        vesselnessType->Update();
        OutputImageType::Pointer outObject = vesselnessType->GetOutput();
        ITKImageDouble outDouble(outObject);
        ITKImageFloat outFloat(outDouble.domain());
        for (const Point &p : outDouble.domain()) {
            outFloat.setValue(p, (float) outDouble(p));
        }
        string outRadiusName = outname.substr(0, lastindex) + "_" + to_string(i) + extension;
        ITKWriter<ITKImageFloat>::exportITK(outRadiusName, outFloat);
    }
    DGtal::trace.endBlock();

    viewer << Viewer3D<>::updateDisplay;
    app.exec();

    return 0;
}
