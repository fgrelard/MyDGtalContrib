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

#include "geometry/CurveProcessor.h"
#include "geometry/CurveDecomposition.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template <typename Container>
void exportToSDP(std::string outputFileName, const Container& points) {
    typedef typename Container::value_type Point;
    ofstream fout;
    fout.open(outputFileName.c_str());
    for (const Point& p : points) {
        for (int i = 0; i < Point::dimension; i++)
            fout << p[i] << " ";
        fout <<std::endl;
    }
    fout.close();
}



int main( int  argc, char**  argv )
{
    typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (curve)")
        ("output,o", po::value<std::string>(), "outputfilename (generate _edge and _vertex files")
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



    Image curve = VolReader<Image>::importVol(curveFilename);
    Z3i::Domain domainCurve = curve.domain();
    Z3i::DigitalSet setCurve(domainCurve);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setCurve, curve,
                                                  thresholdMin-1, thresholdMax);

    CurveProcessor<Z3i::DigitalSet> curveProcessor(setCurve);
    auto branchingPoints = curveProcessor.branchingPoints();
    CurveDecomposition<Z3i::DigitalSet> curveDecomposition(setCurve, branchingPoints);
    auto graph = curveDecomposition.branchDecomposition();
    Z3i::DigitalSet vertex(setCurve.domain());
    for (const Z3i::DigitalSet& branch : graph) {
        CurveProcessor<Z3i::DigitalSet> curveProc(branch);
        Z3i::DigitalSet endPoints = curveProc.endPoints();
        vertex.insert(*endPoints.begin());
        vertex.insert(*std::next(endPoints.begin(), endPoints.size() - 1));
    }
    std::vector<Z2i::Point> edges;
    for (const Z3i::DigitalSet& branch : graph) {
        CurveProcessor<Z3i::DigitalSet> curveProc(branch);
        Z3i::DigitalSet endPoints = curveProc.endPoints();
        Z3i::Point first = *endPoints.begin();
        Z3i::Point end = *std::next(endPoints.begin(), endPoints.size() - 1);
        unsigned int posFirst = std::distance(vertex.begin(), vertex.find(first));
        unsigned int posEnd = std::distance(vertex.begin(), vertex.find(end));
        Z2i::Point v(posFirst, posEnd);
        edges.push_back(v);
    }
    exportToSDP(outputFilename + "_vertex.sdp", vertex);
    exportToSDP(outputFilename + "_edges.sdp", edges);
    return 0;


}
