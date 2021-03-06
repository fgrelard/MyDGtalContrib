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
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "geometry/PointUtil.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


Z3i::DigitalSet readSkelab(string filename)
{

	FILE *f = fopen(filename.c_str(),"r");
	vector<Z3i::Point> points;
	if (f)
	{
		fscanf(f, "ID Cx Cy Cz RADIUS #NEIGHBORS NEIGHBORS_LIST\n");

		int count;
		fscanf(f,"%d\n",&count);
		trace.info() << count << endl;
		for (int i = 0; i < count; ++i)
		{

			int        id;        // unique ID
			float      radius;    // maximal ball radius
			int        valence;   // number of adjcent points
			double a,b,c;
			fscanf(f, "%d %lf %lf %lf %f %d ", &id, &a, &b, &c, &radius, &valence);
//                    std::cout << id << ") " << coord[0] << ", "<< coord[1] << ", "<< coord[2] << endl;
			Z3i::Point coord((int)a, (int)b, (int)c);
			points.push_back(coord);
			for (int j = 0; j < valence; ++j)
			{
				int neigh;

				fscanf(f, "%d ", &neigh);

			}
		}
		fclose(f);

		Z3i::Domain boundingBox = PointUtil::computeBoundingBox<Z3i::Domain>(points);
		trace.info() << boundingBox << endl;
		Z3i::DigitalSet aSet(boundingBox);
		aSet.insert(points.begin(), points.end());
		return aSet;
	}
	return Z3i::DigitalSet(Z3i::Domain(Z3i::Point(0,0,0), Z3i::Point(0,0,0)));
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
	Z3i::Point translationVector(0, 0, 0);
	Metric l2;

	trace.info() << "read file" << endl;
	Z3i::DigitalSet aSet = readSkelab(inputFilename);
	trace.info() << aSet.size() << endl;
	Image image = ImageFromSet<Image>::create(aSet, 1);
	GenericWriter<Image>::exportFile(outputFilename, image);
	return 0;
}
