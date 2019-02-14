#include <iostream>
#include <QtWidgets/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/writers/ITKWriter.h"
#include "geometry/CurveProcessor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
        typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;

        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("input,i", po::value<std::string>(), "vol file (vol)")
                ("output,o", po::value<std::string>(), "output ensured connectivity")
                ("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
                ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
                ;

        bool parseOK=true;
        po::variables_map vm;
        try{
                po::store(po::parse_command_line(argc, argv, general_opt), vm);
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

        string volFilename = vm["input"].as<std::string>();
        string outputFilename = vm["output"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();



        Image vol = VolReader<Image>::importVol(volFilename);
        Z3i::Domain domainVol = vol.domain();
        Z3i::DigitalSet setVol(domainVol);
        SetFromImage<Z3i::DigitalSet>::append<Image> (setVol, vol,
                                                      thresholdMin-1, thresholdMax);

        Z3i::Object26_6 obj(Z3i::dt26_6, setVol);

        Z3i::DigitalSet setVolThinned(setVol);

        for (const Z3i::Point& p : domainVol) {
            if (vol(p) < thresholdMin) continue;
            if (obj.isSimple(p) && obj.neighborhoodSize(p) > 1) {
                setVolThinned.erase(p);
            }
            obj = Z3i::Object26_6(Z3i::dt26_6, setVolThinned);
        }

        Image outImage(setVol.domain());
        DGtal::imageFromRangeAndValue(setVolThinned.begin(), setVolThinned.end(), outImage, 10);
        VolWriter<Image>::exportVol(outputFilename, outImage);
        return 0;


}
