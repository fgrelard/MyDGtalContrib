#include <iostream>

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
#include "vcm/OrthogonalPlaneEstimator.h"
#include "shapes/DigitalPlane.h"
#include "geometry/DigitalPlaneProcessor.h"
#include "geometry/CurveProcessor.h"
#include "geometry/CurveDecomposition.h"
using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
        typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
        typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> OrthoPlaneEstimator;
        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("curve,c", po::value<std::string>(), "vol file (curve)")
                ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
                ("thresholdMin,m", po::value<int>()->default_value(1), "minimum threshold for binarization")
                ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
                ("radiusInside,R", po::value<double>()->default_value(7), "radius of the ball inside voronoi cell")
                ("radiusNeighbour,r", po::value<double>()->default_value(5), "radius of the ball for the neighbourhood")
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
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();
        double R = vm["radiusInside"].as<double>();
        double r = vm["radiusNeighbour"].as<double>();

        Image volume = VolReader<Image>::importVol(inputFilename);
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
        Z3i::DigitalSet endp = curveProcessor.endPoints();
        Z3i::DigitalSet branch = curveProcessor.branchingPoints();
        CurveDecomposition<Z3i::DigitalSet> curveDecompo(setCurve, branch);
        std::vector<Z3i::Point> orderedCurve;
        if (endp.size() > 0)
            orderedCurve = curveDecompo.curveTraversalForGraphDecomposition(*endp.begin());
        else
            orderedCurve = curveDecompo.curveTraversalForGraphDecomposition(*setCurve.begin());


        QApplication application(argc,argv);
        Viewer3D<> viewer;
        viewer.show();


        int sliceNumber = 0;
        int moduloFactor = 10;
        // Compute VCM and diagonalize it.
        viewer.setFillColor(Color::Gray);
        viewer.setFillTransparency(255);

        Z3i::Point current(15,15,54); //it->getPoint();
        Z3i::RealVector normal(-.5226, -0.5226, -0.6736);
        Z3i::RealVector normal2(-.5226, -1.2226, -0.2736);
        DigitalPlane<Z3i::Space> plane(current,normal.getNormalized());
        DigitalPlane<Z3i::Space> plane2(current,normal2.getNormalized());
        DigitalPlaneProcessor<Z3i::Space> planeProc(plane);
        std::vector<Z3i::RealVector> points = planeProc.planeToQuadrangle();
        DigitalPlaneProcessor<Z3i::Space> planeProc2(plane2);

        std::vector<Z3i::RealVector> points2 = planeProc2.planeToQuadrangle();
        double f = 25.0;

        viewer.setLineColor(Color::Blue);
        viewer.setFillColor(Color::Blue);
        viewer.setFillTransparency(150);
        viewer.addQuad(current + points[0] * f,
                       current + points[1] * f,
                       current + points[2] * f,
                       current + points[3] * f);

        viewer.setLineColor(Color::Red);
        viewer.setFillColor(Color::Red);
        viewer.addQuad(current + points2[0] * f,
                       current + points2[1] * f,
                       current + points2[2] * f,
                       current + points2[3] * f);

        typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D;
        Image2D  image1 = planeProc.sliceFromPlane(volume, 100);
        Image2D image2 = planeProc2.sliceFromPlane(volume, 100);

        DGtal::GenericWriter<Image2D>::exportFile("/home/florent/trash/plane1.pgm", image1);
        DGtal::GenericWriter<Image2D>::exportFile("/home/florent/trash/plane2.pgm", image2);

        for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
            if (volume(*it) >= thresholdMin)
                viewer << CustomColors3D(Color(0,0,255,20), Color(0,0,255,20))<<*it;
        }


        viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setVolume;
        viewer << Viewer3D<>::updateDisplay;
        application.exec();
        return 0;


}
