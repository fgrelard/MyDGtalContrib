#include "ShapeDescriptor.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

using namespace DGtal;
using namespace std;
void testMajorAxis() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
        std::string home = getenv("HOME");
        std::string path = home + "/test_img/thskeleton_boudin.vol";
        Image image = GenericReader<Image>::import(path);
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        ShapeDescriptor<Z3i::DigitalSet> setProcessor(setVolume);
        pair<Z3i::Point, Z3i::Point> pair = setProcessor.majorAxis();
        trace.info() << pair.first << " " << pair.second << endl;
}


void testLengthMajorAxis() {
        typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;

        std::string home = getenv("HOME");
        std::string path = home + "/test_img/thskeleton_boudin.vol";
        Image image = GenericReader<Image>::import(path);
        Z3i::DigitalSet setVolume(image.domain());
        SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
        ShapeDescriptor<Z3i::DigitalSet> setProcessor(setVolume);
        double radius = setProcessor.lengthMajorAxis();
        trace.info() << radius << endl;
        const Color CURVE3D_COLOR( 100, 100, 140, 128 );
}

void testCenterOfMass() {
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;

	Image image = GenericReader<Image>::import("/home/florent/test_img/slices/boudin/slice_1.pgm");
	ShapeDescriptor<Z2i::DigitalSet> stats(Z2i::DigitalSet(image.domain()));
	Z2i::Point point = stats.centerOfMass(image);
	trace.info() << point << endl;
}

void testCenterOfMass3D() {
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;

	Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");
	Z3i::DigitalSet set(image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image>(set, image, 0, 255);
	ShapeDescriptor<Z3i::DigitalSet> stats(set);
	Z3i::Point point = stats.extractCenterOfMass();
	trace.info() << point << endl;
}

void testComputeCovarianceMatrix() {
	typedef Eigen::MatrixXd MatrixXd;
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;
	typedef DGtal::Z2i::RealPoint RealPoint;

	Image image = GenericReader<Image>::import("/home/florent/test_img/slices/boudin/slice_1.pgm");
	ShapeDescriptor<Z2i::DigitalSet> stats(Z2i::DigitalSet(image.domain()));
	MatrixXd matrixCovariance = stats.computeCovarianceMatrixImage(image);
	RealPoint vectorZero = stats.extractEigenVector(matrixCovariance, 0);
	RealPoint vectorOne = stats.extractEigenVector(matrixCovariance, 1);
	trace.info() << vectorZero << " " << vectorOne << endl;
}

int main() {
    testMajorAxis();
	testCenterOfMass();
	testCenterOfMass3D();
	testComputeCovarianceMatrix();
	return 0;
}
