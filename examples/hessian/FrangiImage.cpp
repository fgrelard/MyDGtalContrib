
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
#include <hessian/DiscreteHessianFunction.h>
#include <DGtal/io/readers/GenericReader.h>
#include <hessian/FrangiVesselness.h>
#include <DGtal/io/writers/ITKWriter.h>


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

struct DoubleToFloatFunctor {
    float operator()(double value) {
        return static_cast<float>(value);
    }
};


int main(int argc, char **argv) {
    typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
    typedef DGtal::ImageContainerByITKImage<Z3i::Domain, unsigned char> DGtalITKImage;
    typedef DiscreteHessianFunction<DGtalITKImage> Hessian;
    typedef typename Hessian::OutputImage HessianImage;
    typedef DGtal::FrangiVesselness<HessianImage> Vesselness;
    typedef typename Vesselness::OutputImage VesselnessImage;
    typedef ImageContainerBySTLVector<Z3i::Domain, float> WritableImage;


    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "vol file (input volume)")
            ("output,o", po::value<std::string>(), "vol file (output volume)");

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

    DGtalITKImage volume = DGtal::ITKReader<DGtalITKImage>::importITK(inputFilename);


    Hessian hessian(volume, 1.0);
    DGtal::trace.beginBlock("Computing hessian");
    Hessian::OutputImage hessianImage = hessian.computeHessian();
    DGtal::trace.endBlock();

    DGtal::trace.beginBlock("Computing frangi");
    Vesselness vesselness(hessianImage);
    VesselnessImage vesselnessImage = vesselness.computeVesselness();
    WritableImage out(vesselnessImage.domain());
    for (const Z3i::Point &p : out.domain()) {
        out.setValue(p, (float) vesselnessImage(p));
    }
    DGtal::ITKWriter<WritableImage>::exportITK(outputFilename, out);
    DGtal::trace.endBlock();


    return 0;


}
