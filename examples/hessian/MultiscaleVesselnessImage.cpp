
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/io/readers/GenericReader.h>
#include <hessian/FrangiVesselness.h>
#include <DGtal/io/writers/ITKWriter.h>
#include <hessian/HessianRecursiveGaussian.h>
#include <hessian/MultiscaleVesselnessImage.h>


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

struct DoubleToFloatFunctor {
    float operator()(double value) {
        return static_cast<float>(value);
    }
};


int main(int argc, char **argv) {
    using namespace Z3i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> Image;
    typedef DGtal::ImageContainerByITKImage<Domain, unsigned char> DGtalITKImage;
    typedef HessianRecursiveGaussian<Image> Hessian;
    typedef typename Hessian::OutputImage HessianImage;
    typedef DGtal::FrangiVesselness<HessianImage> Vesselness;
    typedef typename Vesselness::OutputImage VesselnessImage;
    typedef ImageContainerBySTLVector<Domain, float> WritableImage;


    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (input volume)")
            ("output,o", po::value<std::string>(), "vol file (output volume)")
    ("sigmaMin,m", po::value<float>()->default_value(1.0), "minimum sigma")
            ("sigmaMax,M", po::value<float>()->default_value(5.0), "maximal sigma")
            ("steps,s", po::value<int>()->default_value(5), "number of step")
            ("alpha,a", po::value<float>()->default_value(0.5), "alpha")
            ("beta,b", po::value<float>()->default_value(1), "beta")
            ("gamma,g", po::value<float>()->default_value(10), "gamma");

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
    string outputFilename = vm["output"].as<std::string>();
    float sigmaMin = vm["sigmaMin"].as<float>();
    float sigmaMax = vm["sigmaMax"].as<float>();
    int steps = vm["steps"].as<int>();
    float alpha = vm["alpha"].as<float>();
    float beta = vm["beta"].as<float>();
    float gamma = vm["gamma"].as<float>();
    Image volume = DGtal::ITKReader<Image>::importITK(inputFilename);


    Hessian hessian(volume, 1.0, true);

    DGtal::trace.beginBlock("Computing multiscale frangi");
    Vesselness vesselness;
    vesselness.setAlpha(alpha);
    vesselness.setBeta(beta);
    vesselness.setGamma(gamma);
    MultiscaleVesselness<Hessian, Vesselness> multiVesselness(hessian, vesselness, sigmaMin, sigmaMax, steps);
    VesselnessImage vesselnessImage = multiVesselness.computeMultiscaleVesselness();
    WritableImage out(vesselnessImage.domain());
    for (const Point &p : out.domain()) {
        out.setValue(p, (float) vesselnessImage(p));
    }
    DGtal::ITKWriter<WritableImage>::exportITK(outputFilename, out);
    DGtal::trace.endBlock();


    return 0;


}
