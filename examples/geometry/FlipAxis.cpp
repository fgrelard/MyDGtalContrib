#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "first vol (reference) .vol" )
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
    string outputFilename = vm["output"].as<std::string>();



    typedef Z3i::Space Space;
    typedef Z3i::KSpace KSpace;
    typedef HyperRectDomain<Space> Domain;
    typedef ImageSelector<Domain, unsigned char>::Type Image;
    typedef Z3i::Object26_6 Object;
    Image firstImage = VolReader<Image>::importVol(firstInput);

    Domain firstDomain = firstImage.domain();


    Z3i::DigitalSet setFirst(firstDomain);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setFirst, firstImage,
                                                  0, 255);



    Z3i::DigitalSet flippedSet(firstDomain);
    for (const auto& p : setFirst) {
        Z3i::Point pTranslated = p;
        pTranslated[0] = p[2];
        pTranslated[1] = p[1];
        pTranslated[2] = p[0];
        flippedSet.insert(pTranslated);
    }
    DGtal::trace.info() << flippedSet.size() << std::endl;
    Image outImage(firstImage.domain());

    DGtal::imageFromRangeAndValue(flippedSet.begin(), flippedSet.end(), outImage, 255);
    VolWriter<Image>::exportVol(outputFilename, outImage);
    return 0;
}
