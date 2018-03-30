#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "CannyEdgeDetector.h"
#include "DGtal/io/writers/ITKWriter.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int argc, char **argv )
{
    using namespace DGtal;
    using namespace DGtal::Z3i;

    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayLevelImage;
    typedef CannyEdgeDetector<GrayLevelImage>::OutputImage         OutImage;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
        ("output,o", po::value<std::string>(), "vol file (corresponding volume)")
        ;


    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, (const char *const *) argv, general_opt), vm);
    } catch(const std::exception& ex){
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


    string inputFilename = vm["input"].as<std::string>();
    string outputFilename = vm["output"].as<std::string>();


    GrayLevelImage img  = GenericReader<GrayLevelImage>::import( inputFilename );
    auto domain = img.domain();

    CannyEdgeDetector<GrayLevelImage> ced(img, 1);
    OutImage out = ced.getOutput();

    ITKWriter<OutImage>::exportITK(outputFilename, out);


    return 0;
}
