#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <algorithm>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include <Eigen/Dense>
#include "DGtal/math/Histogram.h"
#include "DGtal/math/Statistic.h"
#include "geometry/Distance.h"
#include "DGtal/images/CImage.h"

template <typename Container>
class Statistics {

public:
	Statistics() = delete;
	Statistics(const Container& aData) : myData(aData) {}

	double mean();
	double stddev();

	template <typename Image>
	typename Container::Space::RealVector centerOfMass(const Image& image);

	typename Container::Space::RealVector extractCenterOfMass();

	typename Container::Space::RealVector computeNormalFromLinearRegression();

    typename Container::Space::RealVector computeNormalFromCovarianceMatrix();

	template <typename Matrix>
	typename Container::Space::RealVector extractEigenVector(const Matrix& m, int colNumber);

	template <typename Matrix>
	typename Container::Space::RealVector extractEigenValue(const Matrix& m, int colNumber);

	template <typename Matrix>
	Matrix computeCovarianceMatrix();

	template <typename Matrix, typename Image2D>
	Matrix computeCovarianceMatrixImage(const Image2D& image);

	double otsuThreshold();

	double unimodalThresholding();

private:
	Container myData;
};

template <typename Container>
double Statistics<Container>::mean() {
	return std::accumulate( myData.begin(), myData.end(), 0.0f )/ myData.size();
}

template <typename Container>
double Statistics<Container>::stddev() {
	std::vector<double> zero_mean( myData );
	transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(), bind2nd( std::minus<double>(), mean(myData) ) );

	double deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
	deviation = sqrt( deviation / ( myData.size() - 1 ) );
	return deviation;
}


template <typename Container>
template <typename Image>
typename Container::Space::RealVector Statistics<Container>::centerOfMass(const Image& image) {
    BOOST_CONCEPT_ASSERT(( DGtal::concepts::CImage< Image > ));
	double m000 = 0;
	std::vector<double> masses(Container::Space::RealVector::dimension, 0.0);

	for (auto it = myData.domain().begin(), ite = myData.domain().end(); it != ite; ++it) {
	    typename Container::Space::Point current = *it;
		m000 += image(current);
		for (typename Container::Space::RealVector::Dimension i = 0; i < Container::Space::RealVector::dimension; i++) {
			masses[i] += current[i] * image(current);
		}
	}
	if (m000 != 0) {
		typename Container::Space::RealVector v;
		for (typename Container::Space::RealVector::Dimension i = 0; i < Container::Space::RealVector::dimension; i++)
			v[i] = masses[i] * 1.0 / m000;
		return v;
	}

	return typename Container::Space::RealVector();
}

template <typename Container>
typename Container::Space::RealVector Statistics<Container>::extractCenterOfMass() {
	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

	if (myData.size() != 0) {
		typedef typename DGtal::ImageSelector<typename Container::Domain, unsigned char>::Type Image;
		Image image = DGtal::ImageFromSet<Image>::create(myData, 150);
		typename Container::RealVector centerOfMass = Statistics::centerOfMass(image);
		return centerOfMass;
	}
	return typename Container::RealVector();
}


/**
 * Computes the normal of a plane from a set of points
 * Method : linear regression
 */
template <typename Container>
typename Container::Space::RealVector Statistics<Container>::computeNormalFromLinearRegression() {
	typedef Eigen::Matrix<double, Eigen::Dynamic, 3> MatrixXi;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
	unsigned int size = myData.size();
	MatrixXi A(size, 3);
	VectorXi b = VectorXi::Zero(size, 1);

	for (int i = 0; i < size; i++) {
		A(i, 0) = (double)myData[i][0]*1.0;
		A(i, 1) = (double)myData[i][1]*1.0;
		A(i, 2) = 1.0;
		b(i, 0) = (double)myData[i][2]*1.0;
	}
	Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
	typename Container::Space::RealVector normal;
	normal[0] = x(0, 0);
	normal[1] = x(1, 0);
	normal[2] = -1;
	return normal.getNormalized();
}

/**
 * Computes the normal of a plane from a set of points
 * Method : covariance matrix
 */
template <typename Container>
typename Container::Space::RealVector Statistics<Container>::computeNormalFromCovarianceMatrix() {
	typedef Eigen::MatrixXd MatrixXd;

	unsigned int size = myData.size();
	if (size < 2) return (Container::Space::RealVector::zero);

	MatrixXd A(size, 3);
	for (int i = 0; i < size; i++) {
		A(i, 0) = (double)myData[i][0] * 1.0;
		A(i, 1) = (double)myData[i][1] * 1.0;
		A(i, 2) = (double)myData[i][2] * 1.0;
	}
	MatrixXd centered = A.rowwise() - A.colwise().mean();
	MatrixXd cov = (centered.adjoint() * centered) / double(A.rows() - 1);
	Eigen::SelfAdjointEigenSolver<MatrixXd> eig(cov);
	typename Container::Space::RealVector normal;
	auto veigen = eig.eigenvectors().col(0);
	normal[0] = veigen[0];
	normal[1] = veigen[1];
	normal[2] = veigen[2];
	return normal;
}

template <typename Container>
template <typename Matrix>
Matrix Statistics<Container>::computeCovarianceMatrix() {
	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

	typedef typename Container::ConstIterator ConstIterator;
	typedef typename Container::Domain Domain;
	typedef typename Domain::Point Point;

	int dimens = Point::dimension;
	int size = myData.size();
	Matrix A(size, dimens);
	if (size < dimens) return Matrix(0, 0);

	int i = 0;
	for (ConstIterator it = myData.begin(), ite = myData.end();
		 it != ite; ++it) {
		Point point = *it;
		for (int j = 0; j < dimens; j++)
			A(i, j) = (double) point[j] * 1.0;
		i++;
	}
	Matrix centered = A.rowwise() - A.colwise().mean();
	Matrix cov = (centered.adjoint() * centered) / double(A.rows() - 1);
    return cov;
}

template <typename Container>
template <typename Matrix, typename Image>
Matrix Statistics<Container>::computeCovarianceMatrixImage(const Image& image) {
	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CImage< Image > ));

	typedef typename Image::Domain Domain;
	typedef typename Domain::Point Point;
	int size = 0;
    Container aSet(image.domain());
	for (typename Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
		 it != ite; ++it) {
		Point point = *it;
		if (image(*it) > 0) {
			size++;
			aSet.insert(*it);
		}
	}
	myData = aSet;
	return computeCovarianceMatrix<Matrix>();
}

template <typename Container>
template <typename Matrix>
typename Container::Space::RealVector Statistics<Container>::extractEigenVector(const Matrix& m, int colNumber) {
	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(m);
	typename Container::Space::RealVector vector;
	auto veigen = eig.eigenvectors().col(colNumber);
	for (typename Container::Space::RealVector::Dimension i = 0; i < Container::Space::RealVector::dimension; i++) {
		vector[i] = veigen[i];
	}
	return vector;
}

template <typename Container>
template <typename Matrix>
typename Container::Space::RealVector Statistics<Container>::extractEigenValue(const Matrix& m, int colNumber) {
	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDigitalSet< Container > ));

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(m);
	typename Container::Space::RealVector vector;
	auto veigen = eig.eigenvalues().col(colNumber);
	for (typename Container::Space::RealVector::Dimension i = 0; i < Container::Space::RealVector::dimension; i++) {
		vector[i] = veigen[i];
	}
	return vector;
}


template <typename Container>
double Statistics<Container>::otsuThreshold() {
	using namespace DGtal;
	double  proba = 0;                // first order cumulative
    double  mu = 0;                // second order cumulative
    double  mean = 0;               // total mean level
	double    threshold = 0;        // optimal threshold value
	double max = 0.0;

	Statistic<double> stats;
	stats.addValues( myData.begin(), myData.end() );
	stats.terminate(); // stats are computed.

	Histogram<double>* hist = new Histogram<double>();
	hist->init( Histogram<double>::SquareRoot, stats );
	hist->addValues( myData.begin(), myData.end() );
	hist->terminate();
	double myWidth = ( stats.max() - stats.min() ) / hist->size() ;
	double myBin = stats.min();
	for (int i=0; i< hist->size(); i++) {
		myBin += myWidth;
//		std::cout << myBin << " " << hist->pdf(i) << endl;
		mean+= ((double) i / hist->size()) * hist->pdf(i);
	}
	for (int i = 0; i < hist->size(); i++) {
		proba += hist->pdf(i);
		mu += ((double)i/hist->size()) * hist->pdf(i);
		double currentValue =  pow((mean * proba - mu), 2) * proba * (1 - proba);
		if (currentValue > max) {
			max = currentValue;
			threshold = ((double)i/hist->size());
		}

	}

	return threshold;
}

template <typename Container>
double Statistics<Container>::unimodalThresholding() {
	using namespace DGtal;
	Statistic<double> stats;
	stats.addValues( myData.begin(), myData.end() );
	stats.terminate(); // stats are computed.

	Histogram<double>* hist = new Histogram<double>();
	hist->init( Histogram<double>::SquareRoot, stats );
	hist->addValues( myData.begin(), myData.end() );
	hist->terminate();
	double myWidth = ( stats.max() - stats.min() ) / hist->size() ;
	Z2i::RealPoint maxPeak(0,0);
	for (int i = 1; i < hist->size(); i++) {
//		cout << i*myWidth+stats.min() << " " << hist->pdf(i) << endl;
		if (hist->pdf(i) > maxPeak[1])
			maxPeak = Z2i::RealPoint(i*myWidth, hist->pdf(i));
	}
	Z2i::RealPoint tail(stats.max(), hist->pdf(hist->size()-1));
	Z2i::RealVector directionLine = (tail - maxPeak).getNormalized();
	double maxDistanceOrthogonal = 0.0;
	double threshold = 0.0;

	//Start from maxPeak (origin)
	int begin = maxPeak[0] / myWidth;
	for (int i = begin+1; i < hist->size(); i++) {
		Z2i::RealPoint currentPoint(i * myWidth, hist->pdf(i));
		Z2i::RealVector v = currentPoint - maxPeak;
		Z2i::RealPoint orthogonalProjection = ((v.dot(directionLine)) / (directionLine.dot(directionLine))) * directionLine;

		//Need to change basis (go back to true origin)
		orthogonalProjection += maxPeak;
		double currentOrthogonalDistance = Distance::euclideanDistance(orthogonalProjection, currentPoint);
		if (currentOrthogonalDistance > maxDistanceOrthogonal) {
			maxDistanceOrthogonal = currentOrthogonalDistance;
			threshold = currentPoint[0];
		}
	}
	threshold = threshold + (threshold - maxPeak[0]);
	return threshold;
}


#endif
