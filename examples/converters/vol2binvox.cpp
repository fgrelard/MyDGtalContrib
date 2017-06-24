#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


template <typename TImage>
void write_binvox(std::ostream& out, const TImage& img)
{

//
// read header
//
    Z3i::Point extent = img.extent();
    int width = extent[0];
    int height = extent[1];
    int depth = extent[2];
    DGtal::trace.info() << width <<  " " << height << " " << depth << std::endl;
    out << "#binvox 1" << std::endl;
    out << "dim " << width << " " << height << " " << depth << std::endl;
    out << "data" << std::endl;
    Z3i::Point lowerBound = img.domain().lowerBound();
    Z3i::Point upperBound = img.domain().upperBound();
    for (int i = lowerBound[0]; i < upperBound[0] + 1; i++) {
        for (int j = lowerBound[1]; j < upperBound[1] + 1; j++) {
            for (int k = lowerBound[2]; k < upperBound[2] + 1; k++) {
                Z3i::Point p(i, j, k);
                unsigned char value = img(p);
                unsigned char count = 1;
                if (value)
                    value = (unsigned char) 1;
                out << value << count;
            }
        }
    }
    // for (auto it = img.domain().begin(), ite = img.domain().end(); it != ite; ++it) {
    //     DGtal::trace.info() << *it << std::endl;
    //     unsigned char value = img(*it);
    //     out << value << (unsigned char)1;
    // }
}

int main(int argc, char** argv) {

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
        ("output,o",  po::value<std::string>(), "output itk file" ) ;

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

    string inputFilename = vm["input"].as<std::string>();
    string outputFilename = vm["output"].as<std::string>();
    typedef Z3i::Space Space;
    typedef Z3i::KSpace KSpace;
    typedef HyperRectDomain<Space> Domain;
    typedef ImageSelector<Domain, unsigned char>::Type Image;
    typedef ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric

    typedef EigenDecomposition<3,double> LinearAlgebraTool;

    QApplication application(argc,argv);
    Viewer3D<> viewer;

    Image image = DGtal::VolReader<Image>::importVol(inputFilename);
    std::ofstream out(outputFilename);
    write_binvox(out, image);
    DGtal::trace.info() << "done" << std::endl;

    return 0;
}
