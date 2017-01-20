#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/writers/ITKWriter.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template<typename TPoint>
inline
bool
export2NOFF(std::ostream & out,
			const Z3i::DigitalSet& setSkeleton
			) throw(DGtal::IOException){
	DGtal::IOException dgtalio;
	Z3i::Object26_6 obj(Z3i::dt26_6, setSkeleton);
	try
    {
		out << "ID Cx Cy Cz RADIUS #NEIGHBORS NEIGHBORS_LIST"<< std::endl;
		out << setSkeleton.size()<< std::endl;
		int i = 0;
		for (const TPoint& current : setSkeleton) {
			out << i << " " << current[0] << " " << current[1] << " " << current[2] << " " << 1;
			vector<TPoint> neighbors;
			back_insert_iterator<vector<TPoint>> inserter(neighbors);
			obj.writeNeighbors(inserter, current);
			out << " " << neighbors.size();
			for (const TPoint& n : neighbors) {
				auto iterator = setSkeleton.find(n);
				if (iterator != setSkeleton.end()) {
					int index = std::distance(setSkeleton.begin(), iterator);
					out << " " << index;
				}
			}
			out << std::endl;
			i++;
		}

    }catch( ... )
    {
		trace.error() << "OFF writer IO error on export " << std::endl;
		throw dgtalio;
    }

	return true;
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

	QApplication application(argc,argv);
    Viewer3D<> viewer;
	Z3i::Point translationVector(0, 0, 0);
	Image image = VolReader<Image>::importVol(inputFilename);
	Z3i::DigitalSet set3d (image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image>(set3d, image, 0,255);

	Metric l2;

	ofstream fichier(outputFilename);
	export2NOFF<Z3i::Point>(fichier, set3d);

	return 0;
}
