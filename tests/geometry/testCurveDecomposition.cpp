#include "geometry/CurveDecomposition.h"
#include "geometry/CurveProcessor.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <QtGui/qapplication.h>

using namespace DGtal;
using namespace std;

void testCurveTraversalAndGraph(int argc, char** argv) {
        trace.beginBlock("Test curve traversal");
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/Pruning/PF/skeletonBroncheConnexity.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);
        auto branchingPoints = curveProcessor.branchingPoints();

        CurveDecomposition<Z3i::DigitalSet> curveDecomposition(setVolume, branchingPoints);
        auto curve = curveDecomposition.curveTraversalForGraphDecomposition(*setVolume.begin());
        std::vector<Edge<Z3i::DigitalSet> > graph = curveDecomposition.constructGraph(curve);

        for (const Edge<Z3i::DigitalSet>& edge : graph) {
                int r = rand() % 256, g = rand() % 256, b = rand() % 256;
                viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << edge;
        }
        viewer << Viewer3D<>::updateDisplay;
        app.exec();
        trace.endBlock();
}


void testCurveHierarchicalDecomposition(int argc, char** argv) {
        trace.beginBlock("Test curve hierarchical decomp");
        QApplication app(argc, argv);
        Viewer3D<> viewer;
        viewer.show();

        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        Image image = GenericReader<Image>::import("/home/florent/test_img/Pruning/PF/skeletonBroncheConnexity.vol");
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);

        CurveProcessor<Z3i::DigitalSet> curveProcessor(setVolume);
        auto branchingPoints = curveProcessor.branchingPoints();
        auto endPoints = curveProcessor.endPoints();

        CurveDecomposition<Z3i::DigitalSet> curveDecomposition(setVolume, branchingPoints);
        auto curve = curveDecomposition.curveTraversalForGraphDecomposition(*setVolume.begin());
        std::vector<Edge<Z3i::DigitalSet> > graph = curveDecomposition.constructGraph(curve);
        std::vector<WeightedEdge<Z3i::DigitalSet>* > hierarchicalGraph = curveDecomposition.hierarchicalDecomposition(graph, endPoints);
        for (WeightedEdge<Z3i::DigitalSet>* edge : hierarchicalGraph) {
                int label = edge->getLabel();
                int r = label * 20 % 256;
                int g = label * 50 % 256;
                int b = label * 80 % 256;
                viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << *edge;
        }

        viewer << Viewer3D<>::updateDisplay;
        app.exec();
        trace.endBlock();
}




int main(int argc, char** argv) {
        srand(time(NULL));
        //testCurveTraversalAndGraph(argc, argv);
        testCurveHierarchicalDecomposition(argc, argv);
        return 0;
}
