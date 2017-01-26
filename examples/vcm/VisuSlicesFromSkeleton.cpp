#include <iostream>
#include <QtGui/qapplication.h>
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

#include "geometry/MedialAxis.h"
#include "vcm/PruningOrthogonalPlanes.h"
#include "geometry/DigitalPlaneProcessor.h"
#include "viewer/ViewerSlice.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
        typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
        typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D;
        typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> OrthoPlaneEstimator;
        typedef DistanceTransformation<Z3i::Space, Z3i::DigitalSet, Z3i::L2Metric> DTL2;
        typedef DigitalPlane<Z3i::Space> Plane;

        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("curve,c", po::value<std::string>(), "vol file (curve)")
                ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
                ("volume,v", po::value<std::string>(), "vol file (corresponding volume")
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

        string curveFilename = vm["curve"].as<std::string>();
        string inputFilename = vm["input"].as<std::string>();
        string volumeFilename  = vm["volume"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();

        Image input = VolReader<Image>::importVol(inputFilename);
        Image volume = VolReader<Image>::importVol(volumeFilename);
        Z3i::Domain domainVolume = volume.domain();
        Z3i::DigitalSet setVolume(domainVolume);
        SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                      thresholdMin-1, thresholdMax);

        Image curve = VolReader<Image>::importVol(curveFilename);
        Z3i::Domain domainCurve = curve.domain();
        Z3i::DigitalSet setCurve(domainCurve);
        SetFromImage<Z3i::DigitalSet>::append<Image> (setCurve, curve,
                                                      thresholdMin-1, thresholdMax);

        CurveProcessor<Z3i::DigitalSet> curveProcessor(setCurve);
        Z3i::DigitalSet endPoints = curveProcessor.endPoints();
        Z3i::DigitalSet branchingPoints = curveProcessor.branchingPoints();
        Z3i::Point point = *min_element(setCurve.begin(), setCurve.end(), [&](const Z3i::Point& p1, const Z3i::Point& p2) {
                        return (p1[2] < p2[2]);
                });
        std::vector<Z3i::Point> curveOrdered = CurveDecomposition<Z3i::DigitalSet>(setCurve, branchingPoints).curveTraversalForGraphDecomposition(point);



        DTL2 dt(setVolume.domain(), setVolume, Z3i::l2Metric);
        double radiusVCM = dt(*max_element(domainVolume.begin(), domainVolume.end(), [&](const Z3i::Point& p1, const Z3i::Point& p2) {
                        return (dt(p1) < dt(p2));
                        })) * 2;
        KernelFunction chi(1.0, radiusVCM);
        OrthoPlaneEstimator orthogonalPlaneEstimator(setVolume, chi, 20, 10);
        orthogonalPlaneEstimator.setRadius(radiusVCM);
        std::vector<Plane> planes;
        for (const Z3i::Point& p : curveOrdered) {
                double radius = dt(p) + 2.0;
                orthogonalPlaneEstimator.setRadius(radius);
                Plane plane = orthogonalPlaneEstimator.convergentPlaneAt(p, setVolume, radiusVCM);
                planes.push_back(plane);
        }

        QApplication application(argc,argv);
        ViewerSlice<Z3i::Space, Z3i::KSpace> viewer(planes, input, 100);
        Color color(210,210,210,20);
        viewer << CustomColors3D(Color::Red, Color::Red) << setCurve;
        viewer << CustomColors3D(color, color) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
        viewer.show();
        application.exec();
        return 0;


}
