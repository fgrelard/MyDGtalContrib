#include "Statistics.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"

using namespace DGtal;
using namespace std;

void testCenterOfMass() {
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;

	Image image = GenericReader<Image>::import("/home/florent/test_img/slices/boudin/slice_1.pgm");
	Statistics<Z2i::DigitalSet> stats(Z2i::DigitalSet(image.domain()));
	Z2i::Point point = stats.centerOfMass(image);
	trace.info() << point << endl;
}

void testComputeCovarianceMatrix() {
	typedef Eigen::MatrixXd MatrixXd;
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;
	typedef DGtal::Z2i::RealPoint RealPoint;

	Image image = GenericReader<Image>::import("/home/florent/test_img/slices/boudin/slice_1.pgm");
	Statistics<Z2i::DigitalSet> stats(Z2i::DigitalSet(image.domain()));
	MatrixXd matrixCovariance = stats.computeCovarianceMatrixImage<MatrixXd>(image);
	RealPoint vectorZero = stats.extractEigenVector(matrixCovariance, 0);
	RealPoint vectorOne = stats.extractEigenVector(matrixCovariance, 1);
	trace.info() << vectorZero << " " << vectorOne << endl;
}

int main() {
	testCenterOfMass();
	testComputeCovarianceMatrix();
	return 0;
}
