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
#include "DGtal/io/writers/ITKWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <itkRGBPixel.h>

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

Z3i::DigitalSet difference(const Z3i::DigitalSet& first, const Z3i::DigitalSet& second) {
    Z3i::DigitalSet diff(first.domain());
    for (const Z3i::Point&  p : first) {
        if (second.find(p) == second.end()) {
            diff.insert(p);
        }
    }
    return diff;
}

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("second,s", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
		("first,f",  po::value<std::string>(), "vol file" )
        ("distance,d",  po::value<std::string>(), "vol file" )
        ("output,o",  po::value<std::string>(), "vol file" );

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

    string observedName = vm["first"].as<std::string>();

    string outfilename = vm["output"].as<std::string>();
    string distanceName = vm["distance"].as<std::string>();
	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
    typedef DGtal::GradientColorMap<unsigned char, CMAP_JET > Jet;
    typedef itk::RGBPixel<unsigned char> Pixel;
    typedef itk::Image<Pixel, 3> ITKImage;
    typedef ImageContainerByITKImage<Z3i::Domain, Pixel> OutImage;
    typedef itk::ImageFileWriter< ITKImage > WriterType;
	Image imageObserved = ITKReader<Image>::importITK(observedName);
	Z3i::DigitalSet aSetObserved(imageObserved.domain());
	SetFromImage<Z3i::DigitalSet>::append(aSetObserved, imageObserved, 1, 255);

    Image imageTheoretical(imageObserved.domain());
    if (vm.count("second")) {
        string theoreticalName = vm["second"].as<std::string>();
        imageTheoretical = ITKReader<Image>::importITK(theoreticalName);
    }
	Z3i::DigitalSet aSetTheoretical(imageTheoretical.domain());
	SetFromImage<Z3i::DigitalSet>::append(aSetTheoretical, imageTheoretical, 1, 255);

    Image imageDistance = ITKReader<Image>::importITK(distanceName);
	Z3i::DigitalSet aSetDistance(imageDistance.domain());
	SetFromImage<Z3i::DigitalSet>::append(aSetDistance, imageDistance, 1, 255);


    Z3i::DigitalSet diff = difference(aSetObserved, aSetTheoretical);
    Pixel white;
    white.Set(255, 255, 255);
    Pixel green;
    green.Set(18, 211, 45);
    Pixel black;
    black.Set(0, 0, 0);
    Pixel blue;
    blue.Set(72, 84, 236);
    ITKImage::Pointer outItk = ITKImage::New();
    itk::Size<3> size;
    ITKImage::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    Z3i::RealVector extent = imageObserved.extent();
    DGtal::trace.info() << diff.size() << std::endl;

    size[0] = extent[0];
    size[1] = extent[1];
    size[2] = extent[2];
    ITKImage::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);
    outItk->SetRegions(region);
    outItk->Allocate();
    for (const Z3i::Point& p : aSetObserved.domain()) {
        ITKImage::IndexType index;
        for (int i = 0; i < 3; i++)
            index[i] = p[i];
        if (aSetDistance.find(p) != aSetDistance.end())
            outItk->SetPixel(index, blue);
        else if (diff.find(p) != diff.end())
            outItk->SetPixel(index, green);
        else if (aSetTheoretical.find(p) != aSetTheoretical.end())
            outItk->SetPixel(index, white);
        else
            outItk->SetPixel(index, black);
    }
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outfilename );
    writer->SetInput(outItk);
    writer->Update();
    DGtal::trace.info() << "done " << std::endl;

	return 0;
}
