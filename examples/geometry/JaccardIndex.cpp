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
using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("theoretical,t", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
        ("thresholdMin,m", po::value<double>()->default_value(1), "threshold observed")
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
    double thresholdMin = vm["thresholdMin"].as<double>();

	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;

	Image imageObserved = ITKReader<Image>::importITK(observedName);
	Z3i::DigitalSet aSetObserved(imageObserved.domain());
	SetFromImage<Z3i::DigitalSet>::append(aSetObserved, imageObserved, thresholdMin, 255);

	Image imageTheoretical = ITKReader<Image>::importITK(theoreticalName);
	Z3i::DigitalSet aSetTheoretical(imageTheoretical.domain());
	SetFromImage<Z3i::DigitalSet>::append(aSetTheoretical, imageTheoretical, 1, 255);


    DGtal::trace.info() << aSetObserved.size() << " " << aSetTheoretical.size() << std::endl;
	Z3i::DigitalSet intersection(imageTheoretical.domain());
	for (const Z3i::Point& t : aSetTheoretical) {
		for( const Z3i::Point& o : aSetObserved) {
			if (t == o)
				intersection.insert(t);
		}
	}


    DGtal::trace.info() << intersection.size() << std::endl;
    Z3i::DigitalSet unionSet(aSetObserved);
    unionSet.insert(aSetTheoretical.begin(), aSetTheoretical.end());

	double precision = (1.0 * intersection.size() / aSetObserved.size());
	double recall = (1.0 * intersection.size() / aSetTheoretical.size());

	double fmeasure = (2 * precision * recall) / (precision + recall);

	double distanceMax = 0;
	// for (auto it = aSetObserved.begin(), ite = aSetObserved.end(); it != ite; ++it) {
	// 	Z3i::Point closestPointInTheoretical =  *min_element(aSetTheoretical.begin(), aSetTheoretical.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
	// 				return Z3i::l2Metric(one, *it) < Z3i::l2Metric(two, *it);
	// 			});
	// 	double distance = Z3i::l2Metric(closestPointInTheoretical, *it);
	// 	if (distance > distanceMax)
	// 		distanceMax = distance;
	// }
    double jaccard = intersection.size() * 1.0 / unionSet.size();
    double diff = aSetObserved.size()*1.0 - aSetTheoretical.size();
//	trace.info() << "h-distance " <<  distanceMax << endl;
	trace.info() << "f-measure " << fmeasure << endl;
    trace.info() << "jaccard " << jaccard << std::endl;
    trace.info() << "precision " << precision << std::endl;
    trace.info() << "recall " << recall << std::endl;
    trace.info() << "diff " << diff << std::endl;
	return 0;
}
