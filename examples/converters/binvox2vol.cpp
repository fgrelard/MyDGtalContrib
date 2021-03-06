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


typedef unsigned char byte;



static int version;
static int depth, height, width;
static int size;
static byte *voxels = 0;
static float tx, ty, tz;
static float scale;



Z3i::DigitalSet read_binvox(string filespec, const Z3i::DigitalSet& set)
{
	Z3i::DigitalSet voxelSet(set.domain());	
	ifstream *input = new ifstream(filespec.c_str(), ios::in | ios::binary);

//
// read header
//
	string line;
	*input >> line;  // #binvox
	if (line.compare("#binvox") != 0) {
		cout << "Error: first line reads [" << line << "] instead of [#binvox]" << endl;
		delete input;
	}
	*input >> version;
	cout << "reading binvox version " << version << endl;

	depth = -1;
	int done = 0;
	while(input->good() && !done) {
		*input >> line;
		if (line.compare("data") == 0) done = 1;
		else if (line.compare("dim") == 0) {
			*input >> depth >> height >> width;
		}
		else if (line.compare("translate") == 0) {
			*input >> tx >> ty >> tz;
		}
		else if (line.compare("scale") == 0) {
			*input >> scale;
		}
		else {
			cout << "  unrecognized keyword [" << line << "], skipping" << endl;
			char c;
			do {  // skip until end of line
				c = input->get();
			} while(input->good() && (c != '\n'));

		}
	}
	if (!done) {
		cout << "  error reading header" << endl;
	}
	if (depth == -1) {
		cout << "  missing dimensions in header" << endl;
	}

	size = width * height * depth;
	voxels = new byte[size];
	if (!voxels) {
		cout << "  error allocating memory" << endl;
	}

//
// read voxel data
//
	byte value;
	byte count;
	int index = 0;
	int end_index = 0;
	int nr_voxels = 0;
  
	input->unsetf(ios::skipws);  // need to read every byte now (!)
	*input >> value;  // read the linefeed char
	
	while((end_index < size) && input->good()) {
		*input >> value >> count;
        if (input->good()) {
			end_index = index + count;
			if (end_index > size) break;
			for(int i=index; i < end_index; i++)
			{
				voxels[i] = value;
				if (value) {
                    int x = i / (width*height);
                    int y = (i / width) % height;
                    int z = i % width;
					voxelSet.insert(Z3i::Point(x,y,z));
				}
			}
      
			if (value) {
				nr_voxels += count;

			}
			index = end_index;
		}  // if file still ok
    
	}  // while

	input->close();
	cout << "  read " << nr_voxels << " voxels" << endl;

	return voxelSet;

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

	Z3i::DigitalSet set3d (Domain(Z3i::Point(0,0,0), Z3i::Point(512,512,512)));
	Z3i::DigitalSet voxelSet = read_binvox(inputFilename, set3d);
	trace.info() << voxelSet.size() << endl;
	Image out = ImageFromSet<Image>::create(voxelSet, 255);
	VolWriter<Image>::exportVol(outputFilename, out);


	return 0;
}
