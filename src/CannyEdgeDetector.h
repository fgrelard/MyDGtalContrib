#include <cmath>
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/math/Statistic.h"

template <typename Image>
class CannyEdgeDetector {
public:
    typedef typename Image::Domain Domain;
    typedef typename Domain::Space Space;
    typedef typename Space::Point Point;
    typedef typename Space::RealVector RealVector;
    typedef DGtal::ImageContainerBySTLVector<Domain, RealVector> DirImage;
    typedef DGtal::ImageContainerBySTLVector<Domain, float> OutputImage;
    typedef DGtal::MetricAdjacency<Space, 3> MetricAdjacency;

public:
    CannyEdgeDetector() = delete;
    CannyEdgeDetector(const Image& aImage, double aSigma = 0.0) : myImage(aImage), mySigma(aSigma), myDirectionGradient(aImage.domain()) {}
    CannyEdgeDetector(const CannyEdgeDetector& ced) : myImage(ced.myImage), mySigma(ced.mySigma), myDirectionGradient(myImage.domain()){}

public:
    OutputImage getOutput();

private:
    OutputImage smoothedImage();
    float gaussianBlur(int x);

    float gradientValue(const OutputImage& image, const Point& p);
    bool onBorder(const Point& p);

    OutputImage suppressNonMaximumGradient(const OutputImage& out);

private:
    Image myImage;
    double mySigma;


    DirImage myDirectionGradient;
};



template <typename Image>
typename CannyEdgeDetector<Image>::OutputImage
CannyEdgeDetector<Image>::getOutput() {
    Domain domain = myImage.domain();
    OutputImage out(domain);
    OutputImage smoothed = smoothedImage();
    DGtal::trace.info() << "Smoothed" << std::endl;
    for (const Point& p : domain) {
        if (onBorder(p)) {
            out.setValue(p, 0);
        }
        else {
            float value = gradientValue(smoothed, p);
            out.setValue(p, value);
        }
    }
    DGtal::trace.info() << "Done." << std::endl;
    return suppressNonMaximumGradient(out);
}


template <typename Image>
typename CannyEdgeDetector<Image>::OutputImage
CannyEdgeDetector<Image>::smoothedImage() {
    Domain domain = myImage.domain();
    OutputImage out(domain);
    for( const Point& p : domain) {
        out.setValue(p, myImage(p));
    }
    for (int i = 0; i < Point::dimension; i++) {
        for (const Point& p : domain) {
            float smoothed = 0;
            for (int j = -2; j <= 2; j++) {
                Point pn = p;
                pn[i] += j;
                if (!domain.isInside(pn)) {
                    smoothed = 0;
                    break;
                }
                smoothed += gaussianBlur(j) * out(pn);
            }
            out.setValue(p, smoothed);
        }
    }
    return out;
}



template <typename Image>
float CannyEdgeDetector<Image>::gaussianBlur(int x) {
    return (1.0f / (std::sqrt(2 * M_PI) * mySigma)) * std::exp(-(x*x)/(2*mySigma*mySigma));
}

template <typename Image>
bool CannyEdgeDetector<Image>::onBorder(const Point& p) {
    Point up = myImage.domain().upperBound();
    Point low = myImage.domain().lowerBound();

    return (p[0] == up[0] || p[0] == low[0] ||
            p[1] == up[1] || p[1] == low[1] ||
            p[2] == up[2] || p[2] == low[2]);
}

template <typename Image>
float CannyEdgeDetector<Image>::gradientValue(const OutputImage& image, const Point &p) {
    double g2 = 0;

    double greatestMag = 0;
    int index = 0;
    RealVector dir;
    for (int i = 0; i < Point::dimension; i++) {
        Point pp = p;
        Point pn = p;
        pp[i] -= 1;
        pn[i] += 1;
        double gi = -image(pp) + image(pn);

        g2 += gi * gi;
        dir[i] = gi;
        // if (std::abs(gi) > greatestMag) {
        //     greatestMag = std::abs(gi);
        //     index = i;
        // }
    }
    myDirectionGradient.setValue(p, dir.getNormalized());
    return std::sqrt(g2);
}

template <typename Image>
typename CannyEdgeDetector<Image>::OutputImage
CannyEdgeDetector<Image>::suppressNonMaximumGradient(const OutputImage& out) {
    Domain domain = out.domain();
    MetricAdjacency adj;
    OutputImage outNonMax(domain);

    for (const Point& p : domain) {
        outNonMax.setValue(p, out(p));
    }
    for (const Point& p : domain) {

        RealVector dir = myDirectionGradient(p);
        Point pp = p - dir*2;
        Point pn = p + dir*2;
        // pp[dim] -= 1;
        // pn[dim] += 1;
        if (domain.isInside(pn) && domain.isInside(pp) && (out(pn) > out(p) || out(pp) > out(p))) {
            DGtal::Statistic<double> previousValues;
            DGtal::Statistic<double> nextValues;
            for (int i = 2; i < 7; i++) {
                Point pp = p - dir*i;
                Point pn = p + dir*i;
                if (domain.isInside(pp))
                    previousValues.addValue(myImage(pp));
                if (domain.isInside(pn))
                    nextValues.addValue(myImage(pn));
            }
            previousValues.terminate();
            nextValues.terminate();
            // if (previousValues.mean() < nextValues.mean())
            //     outNonMax.setValue(p,100);
            // else
                outNonMax.setValue(p, 0);
        }

    }
    return outNonMax;
}
