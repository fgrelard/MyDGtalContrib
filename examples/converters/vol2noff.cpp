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
#include "vcm/OrthogonalPlaneEstimator.h"
#include "DGtal/math/linalg/EigenDecomposition.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template<typename TPoint, typename PlaneEstimator>
inline
bool
export2NOFF(std::ostream & out,
			const Z3i::DigitalSet& setSurface,
			const PlaneEstimator& planeEstimator) throw(DGtal::IOException){
	DGtal::IOException dgtalio;

	try
    {
		out << "NOFF"<< std::endl;
		out << "# generated from MeshWriter from the DGTal library"<< std::endl;
		out << setSurface.size()  << " " << 0 << " " << std::endl;
		TPoint lower = setSurface.domain().lowerBound();
		TPoint upper = setSurface.domain().upperBound();
		TPoint middle = (lower + middle) / 2;
		int max_x = upper[0] - lower[0];
		int max_y = upper[1] - lower[1];
		int max_z = upper[2] - lower[2];

		for(auto it = setSurface.begin(), ite = setSurface.end(); it != ite; ++it){
			TPoint p = *it;
			double coord_x = (((double)p[0] - lower[0]) * 2.0 / (upper[0] - lower[0])) -1;
			double coord_y = (((double)p[1] - lower[1]) * 2.0 / (upper[1] - lower[1])) -1;
			double coord_z = (((double)p[2] - lower[2]) * 2.0 / (upper[2] - lower[2])) -1;
			auto plane = planeEstimator.planeAt(p);
			Z3i::RealPoint normal = plane.getPlaneEquation().normal();
			out << coord_x << " " << coord_y << " "<< coord_z << " " << normal[0] << " " << normal[1] << " " << normal[2] <<std::endl;
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
	typedef VoronoiCovarianceMeasure<Space,Metric> VCM;
	typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
	typedef EigenDecomposition<3,double> LinearAlgebraTool;

	QApplication application(argc,argv);
    Viewer3D<> viewer;
	Z3i::Point translationVector(0, 0, 0);
	Image image = VolReader<Image>::importVol(inputFilename);
	Z3i::DigitalSet set3d (image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image>(set3d, image, 0,255);

	Metric l2;
	KernelFunction chi( 1.0, 4 );
	OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> orthoEsti(set3d, chi, 10, 4);

	ofstream fichier(outputFilename);
	export2NOFF<Z3i::Point>(fichier, set3d, orthoEsti);

	return 0;
}
