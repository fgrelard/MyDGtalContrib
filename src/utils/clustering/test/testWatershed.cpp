#include "DGtal/base/Common.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include <QtGui/qapplication.h>
#include "../Watershed.h"

using namespace DGtal;
using namespace std;


void determineParent(const map<Z3i::Point, Watershed<Z3i::Point>::WatershedInformation>& watershedResult,
                     const map<Z3i::Point, double> pointToValue,
                     const Z3i::DigitalSet& setVolume, Viewer3D<>& viewer) {
    Z3i::DigitalSet watershedSet(setVolume.domain());
    Z3i::Object26_6 obj(Z3i::dt26_6, setVolume);
    for (const auto & pairWatershed : watershedResult) {
        Z3i::Point current = pairWatershed.first;
        double label = pairWatershed.second.getLabel();
        if (label == WATERSHED) {
            watershedSet.insert(current);
            viewer << CustomColors3D(Color::Silver, Color::Silver) << current;
        }
    }

    Z3i::Object26_6 watershedObj(Z3i::dt26_6, watershedSet);
    vector<Z3i::Object26_6> cc;
    back_insert_iterator<vector<Z3i::Object26_6> > inserter(cc);
    unsigned int nbCC = watershedObj.writeComponents(inserter);
    for (const Z3i::Object26_6& object : cc) {
        Z3i::DigitalSet surroundingVoxels(object.pointSet().domain());
        for (const Z3i::Point& p : object) {
            vector<Z3i::Point> neighbors;
            back_insert_iterator<vector<Z3i::Point> > inserter(neighbors);
            obj.writeNeighbors(inserter, p);
            surroundingVoxels.insert(neighbors.begin(), neighbors.end());
        }
        map<int, set<Z3i::Point> > decompositionMap;
        for (const Z3i::Point& p : surroundingVoxels) {
            int label = watershedResult.at(p).getLabel();
            if (label == WATERSHED) continue;
            decompositionMap[label].insert(p);
        }

        double maxSize = 0;
        set<Z3i::Point> parent;
        for (const pair<int, set<Z3i::Point> >& pairTube : decompositionMap) {
            double sum = 0;
            set<Z3i::Point> currentSet = pairTube.second;
            double currentSize = currentSet.size();
            if (currentSize > maxSize) {
                maxSize = currentSize;
                parent = currentSet;
            }
        }
        for (const Z3i::Point & p : parent) {
            viewer << CustomColors3D(Color::Red, Color::Red) << p;
        }
    }
}

int testVisu(int argc, char** argv) {
     typedef ImageSelector<Z3i::Domain, double>::Type Image;
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    Image image = ITKReader<Image>::importITK("/home/florent/test_img/Eigenvector/airway.mhd");
    map<Z3i::Point, double> pointToValue;
    for (auto it = image.domain().begin(), ite = image.domain().end(); it != ite; ++it) {
        Z3i::Point p = *it;
        pointToValue[p] = image(p);
    }
    double minVal = min_element(pointToValue.begin(), pointToValue.end(), [&](const pair<Z3i::Point, double>& one,
                                                                              const pair<Z3i::Point, double>& two) {
                                    return (one.second < two.second);
                                })->second;
    trace.info() << minVal << std::endl;
    map<Z3i::Point, double> pointToValueObject;

    for (const auto& pair : pointToValue) {
        Z3i::Point p = pair.first;
        double value = pair.second;
        if (value >= minVal+1) {
            pointToValueObject[p] = value;
        }
    }
    double maxVal = max_element(pointToValueObject.begin(), pointToValueObject.end(), [&](const pair<Z3i::Point, double>& one,
                                                                                          const pair<Z3i::Point, double>& two) {
                                    return (one.second < two.second);
                                })->second;
    trace.info() << maxVal << endl;
    GradientColorMap<double, CMAP_JET > hueShade(minVal+1, maxVal);
    for (const auto& pair : pointToValueObject) {
        viewer << CustomColors3D(hueShade(pair.second), hueShade(pair.second)) << pair.first;
    }
    viewer << Viewer3D<>::updateDisplay;
    return app.exec();
}

int testWatershed(int argc, char** argv) {
    typedef ImageSelector<Z3i::Domain, double>::Type Image;

    Image image = ITKReader<Image>::importITK("/home/florent/test_img/Eigenvector/bronche.mhd");
    map<Z3i::Point, double> pointToValue;
    for (auto it = image.domain().begin(), ite = image.domain().end(); it != ite; ++it) {
        Z3i::Point p = *it;
        pointToValue[p] = image(p);
    }
    double minVal = min_element(pointToValue.begin(), pointToValue.end(), [&](const pair<Z3i::Point, double>& one,
                                                                                        const pair<Z3i::Point, double>& two) {
                                    return (one.second < two.second);
                                })->second;
    trace.info() << minVal << std::endl;
    map<Z3i::Point, double> pointToValueObject;
     Z3i::DigitalSet setVolume(image.domain());
    for (const auto& pair : pointToValue) {
        Z3i::Point p = pair.first;
        double value = pair.second;
        if (value >= minVal+1) {
            pointToValueObject[p] = value;
            setVolume.insert(p);
        }

    }

    double threshold = 0.1;
//    threshold = 0;
    Watershed<Z3i::Point> watershed(pointToValueObject, threshold);
    watershed.compute();
    auto resultWatershed = watershed.getWatershed();
    int bins = watershed.getBins();
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    determineParent(resultWatershed, pointToValueObject, setVolume, viewer);
    const Color CURVE3D_COLOR( 100, 100, 140, 20 );
    // GradientColorMap<int, CMAP_JET > hueShade(0, bins);
    // trace.info() << "Bins= " << bins << endl;
    // for (const auto& complexPoint : resultWatershed) {
    //     viewer << CustomColors3D(hueShade(complexPoint.second.getLabel()), hueShade(complexPoint.second.getLabel())) << complexPoint.first;
    //  }
    viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << setVolume;
    viewer << Viewer3D<>::updateDisplay;
    return app.exec();
}

int main(int argc, char** argv) {
    return testWatershed(argc, argv);
}
