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

Z3i::RealVector naiveBinormal(const Z3i::RealVector& t1,
                              const Z3i::RealVector& t2) {
    if (std::abs(t1.dot(t2)) == 1)
        return Z3i::RealVector::zero;
    else {
        return (t1.crossProduct(t2 - t1).getNormalized());
    }
}

Z3i::RealVector closestVector(const Z3i::RealVector& reference,
                              const Z3i::RealVector& first,
                              const Z3i::RealVector& second) {
    return ((first.cosineSimilarity(reference) < second.cosineSimilarity(reference)) ? first : second);
}

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
        std::vector<Z3i::RealVector> tangents;
        std::vector<Z3i::RealVector> vectors2;
        std::vector<Z3i::RealVector> vectors3;
        for (const Z3i::Point& p : orderedCurve) {
            vcm_r = vcm.measure(chi, p);
            if (vcm_r == nil) continue;
            LinearAlgebraTool::getEigenDecomposition(vcm_r, evec, eval);
            tangents.push_back(evec.column(0));
            vectors2.push_back(evec.column(1));
            vectors3.push_back(evec.column(2));
                        // pToValues[p] = eval[1];
            // double ratio = eval[1] / eval[2];
            // Z3i::RealVector v = evec.column(2);
            // Z3i::RealVector dir  = -p;
            // Z3i::RealVector vector;
            // if (v[2] < 0)
            //     vector = -v;
            // else
            //     vector = v;
            // viewer << CustomColors3D(Color::Red, Color::Red);
            // viewer.addLine(p, p+vector*5);
            // DGtal::trace.info() << eval[2] << std::endl;
        }
        DGtal::trace.endBlock();
        for (int i  = 0; i < tangents.size() - 1; i++) {
            Z3i::Point current = orderedCurve[i];
            Z3i::RealVector t1 = tangents[i];
            Z3i::RealVector t2 = tangents[i+1];
            Z3i::RealVector v2 = vectors2[i];
            Z3i::RealVector v3 = vectors3[i];
            Z3i::RealVector naiveBN = naiveBinormal(t1, t2);
            if (naiveBN == Z3i::RealVector::zero) continue;
            Z3i::RealVector realBN = closestVector(naiveBN, v2, v3);
            viewer << CustomColors3D(Color::Red, Color::Red);
            viewer.addLine(current, current+realBN*5);
        }




        viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setCurve;
        viewer << Viewer3D<>::updateDisplay;
        application.exec();
        return 0;


}
