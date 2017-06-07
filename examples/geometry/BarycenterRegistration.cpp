#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "ShapeDescriptor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,f", po::value<std::string>(), "first vol (reference) .vol" )
        ("input2,s", po::value<std::string>(), "second vol .tif" )
        ("toTranslate,t", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
        ("output,o",  po::value<std::string>(), "output itk file" );


    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        parseOK=false;
        trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
    po::notify(vm);
    if( !parseOK || vm.count("help")||argc<=1)
    {
        std::cout << "Usage: " << argv[0] << " [input]\n"
                  << "Display volume file as a voxel set by using QGLviewer"<< endl
                  << general_opt << "\n";
        return 0;
    }
    if(!vm.count("input"))
    {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }
    string firstInput = vm["input"].as<std::string>();
    string secondInput = vm["input2"].as<std::string>();
    string registerInput = vm["toTranslate"].as<std::string>();
    string outputFilename = vm["output"].as<std::string>();



    typedef Z3i::Space Space;
    typedef Z3i::KSpace KSpace;
    typedef HyperRectDomain<Space> Domain;
    typedef ImageSelector<Domain, unsigned char>::Type Image;
    typedef Z3i::Object26_6 Object;
    Image firstImage = VolReader<Image>::importVol(firstInput);
    Image secondImage = ITKReader<Image>::importITK(secondInput);
    Image toRegister = VolReader<Image>::importVol(registerInput);

    Domain firstDomain = firstImage.domain();
    Domain secondDomain = secondImage.domain();


    Z3i::DigitalSet setFirst(firstDomain);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setFirst, firstImage,
                                                  0, 255);

    Z3i::DigitalSet setSecond(secondDomain);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setSecond, secondImage,
                                                  0, 255);

    Z3i::DigitalSet setToRegister(toRegister.domain());
    SetFromImage<Z3i::DigitalSet>::append<Image> (setToRegister, toRegister,
                                                  0, 255);

    ShapeDescriptor<Z3i::DigitalSet> shape1(setFirst);
    ShapeDescriptor<Z3i::DigitalSet> shape2(setSecond);

    Z3i::Point center1 = shape1.extractCenterOfMass();
    Z3i::Point center2 = shape2.extractCenterOfMass();

    Z3i::Vector vectorTranslation = center1 - center2;

    Z3i::DigitalSet registeredSet(firstDomain);
    for (const auto& p : setToRegister) {
        Z3i::Point pTranslated = p + vectorTranslation;
        registeredSet.insert(pTranslated);
    }
    Image outImage(firstImage.domain());

    DGtal::imageFromRangeAndValue(registeredSet.begin(), registeredSet.end(), outImage, 255);
    VolWriter<Image>::exportVol(outputFilename, outImage);
    return 0;
}
