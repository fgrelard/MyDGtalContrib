#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/math/Statistic.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("theoretical,t", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
		("observed,o",  po::value<std::string>(), "vol file" ) ;

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
	if(!vm.count("theoretical") && !vm.count("observed"))
	{
		trace.error() << " The file name was not defined" << endl;
		return 0;
	}
	string theoreticalName = vm["theoretical"].as<std::string>();
	string observedName = vm["observed"].as<std::string>();


	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, double>::Type Image;

	Image imageObserved = ITKReader<Image>::importITK(observedName);
	Image imageTheoretical = ITKReader<Image>::importITK(theoreticalName);
    DGtal::Statistic<double> stats;
    for (const Z3i::Point& p : imageObserved.domain()) {
        double value = imageObserved(p);
        double valueTh = imageTheoretical(p);

        double diff = std::abs(value - valueTh);
        if (valueTh > 0)
            stats.addValue(diff);


    }

    DGtal::trace.info() << "mean " << stats.mean() << std::endl;
    DGtal::trace.info() << "stddev " << std::sqrt(stats.variance()) << std::endl;
    DGtal::trace.info() << "max " << stats.max() << std::endl;


	return 0;
}
