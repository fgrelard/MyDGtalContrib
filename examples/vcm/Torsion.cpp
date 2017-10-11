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
        typedef VCMAdjustableRadius<Z3i::Space, Z3i::L2Metric> VCM;
        typedef typename VCM::MatrixNN Matrix;
        typedef EigenDecomposition<3, double> LinearAlgebraTool;
        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("curve,c", po::value<std::string>(), "vol file (curve)")
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

        string curveFilename = vm["curve"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();
        double R = vm["radiusInside"].as<double>();
        double r = vm["radiusNeighbour"].as<double>();



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

        double radiusVCM = r;
        KernelFunction chi(1.0, radiusVCM);
        Z3i::L2Metric l2Metric;
        VCM vcm(R, r, l2Metric, true);
        vcm.init(setCurve.begin(), setCurve.end());
        Matrix vcm_r, evec, nil;
        Z3i::RealVector eval;
        DGtal::trace.beginBlock("Matrix image");
        std::map<Z3i::Point, double> pToValues;
        for (const Z3i::Point& p : setCurve) {
            vcm_r = vcm.measure(chi, p);
            if (vcm_r == nil) continue;
            LinearAlgebraTool::getEigenDecomposition(vcm_r, evec, eval);
            pToValues[p] = eval[1];
            double ratio = eval[1] / eval[2];
            Z3i::RealVector v = evec.column(2);
            Z3i::RealVector dir  = -p;
            Z3i::RealVector vector;
            if (v[2] < 0)
                vector = -v;
            else
                vector = v;
            viewer << CustomColors3D(Color::Red, Color::Red);
            viewer.addLine(p, p+vector*5);
            DGtal::trace.info() << eval[2] << std::endl;
        }
        DGtal::trace.endBlock();



        viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setCurve;
        viewer << Viewer3D<>::updateDisplay;
        application.exec();
        return 0;


}
