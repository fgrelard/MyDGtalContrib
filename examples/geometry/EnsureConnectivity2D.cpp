#include <iostream>
#include <QtWidgets/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/ITKWriter.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/CurveProcessor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
        typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image;

        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("input,i", po::value<std::string>(), "vol file (curve)")
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

        string curveFilename = vm["input"].as<std::string>();
        string outputFilename = vm["output"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();



        Image curve = ITKReader<Image>::importITK(curveFilename);
        Z2i::Domain domainCurve = curve.domain();
        Z2i::DigitalSet setCurve(domainCurve);
        SetFromImage<Z2i::DigitalSet>::append<Image> (setCurve, curve,
                                                      thresholdMin-1, thresholdMax);

        Z2i::Object4_8 obj(Z2i::dt4_8, setCurve);
        vector<Z2i::Object4_8> objects;
        back_insert_iterator< vector<Z2i::Object4_8> > inserter(objects);
        obj.writeComponents(inserter);

        double maxSize = 0;
        for (const Z2i::Object4_8& obj : objects) {
            if (obj.size() > maxSize) {
                setCurve = obj.pointSet();
                maxSize = obj.size();
            }
        }


        CurveProcessor<Z2i::DigitalSet> curveProcessor(setCurve);
        Z2i::DigitalSet curveConnectivity = curveProcessor.ensureConnectivityComplementary();

        CurveProcessor<Z2i::DigitalSet> curveProcessor2(curveConnectivity);
        std::vector<Z2i::Point> ordered = curveProcessor2.convertToOrderedCurve();

        for (const auto & p: ordered) {
            DGtal::trace.info() << p[0] << " " << p[1] << std::endl;
        }

        Image outImage(setCurve.domain());
        DGtal::imageFromRangeAndValue(curveConnectivity.begin(), curveConnectivity.end(), outImage, 255);
        ITKWriter<Image>::exportITK(outputFilename, outImage);
        return 0;


}
