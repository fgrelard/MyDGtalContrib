#ifndef __MORPHOMATHS__
#define __MORPHOMATHS__

#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/domains/DomainPredicate.h"

template <typename Image>
class Morphomaths {

public:
	Morphomaths() = delete;
	Morphomaths(const Image& image) : myImage(image) {}
	Morphomaths(const Morphomaths& other) : myImage(other.myImage) {}

public:
	bool process(int x, int y, int nx, int ny);

	Image constructOnePxBorderImage();

	Image erosion(int size);

	Image dilation(int size);

	Image open(int size);

	Image close(int size);

protected:
	Image myImage;

};


template <typename Image>
bool Morphomaths<Image>::process(int x, int y, int nx, int ny) {
	typedef typename Image::Point Point;

	int valueMin = std::numeric_limits<int>::max();
    int valueMax = 0;
	for (int i = -nx; i <= nx; i++) {
		for (int j = -ny; j <= ny; j++) {
			Point current(x+i,y+j);
			if (myImage.domain().isInside(current)) {
				int value = myImage(current);
				if ((i != 0 || j != 0)) {
					if (value < valueMin) valueMin = value;
					if (value > valueMax) valueMax = value;
				}
			}
		}
	}
	return (valueMin != valueMax);
}

template <typename Image>
Image Morphomaths<Image>::constructOnePxBorderImage() {
	typedef typename Image::Domain Domain;
	typedef typename Image::Point Point;

	Domain domain = myImage.domain();

	Image toReturn(Domain(domain.lowerBound() - Point::diagonal(), domain.upperBound() + Point::diagonal()));
	for (auto it = toReturn.domain().begin(), ite = toReturn.domain().end(); it != ite; ++it) {
		if (domain.isInside(*it))
			toReturn.setValue(*it, image(*it));
		else
			toReturn.setValue(*it, 0);
	}

	return toReturn;
}

template <typename Image>
Image Morphomaths<Image>::erosion(int size) {
	typedef typename Image::Domain Domain;
	typedef typename Image::Point Point;

	Domain domain = myImage.domain();
	Image toReturn = constructOnePxBorderImage(myImage);
	myImage = toReturn;
	Point upper = domain.upperBound(), lower = domain.lowerBound();

	int width = upper[0] - lower[0]+1;
	int height = upper[1] - lower[1]+1;

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			bool shouldBeEroded = process(i, j, size, size);
			if (shouldBeEroded) toReturn.setValue(Point(i,j), 0);
		}
	}
	Image out(domain);
	for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
		out.setValue(*it, toReturn(*it));
	}
	return out;
}

template <typename Image>
Image Morphomaths<Image>::dilation(int size) {
	typedef typename Image::Domain Domain;
	typedef typename Image::Point Point;
	typedef DGtal::functors::NotPointPredicate<Image> BackgroundPredicate;

	Domain domain = myImage.domain();
	Image toReturn = constructOnePxBorderImage(myImage);
	myImage = toReturn;
	BackgroundPredicate backgroundPredicate( myImage );


	Point upper = domain.upperBound(), lower = domain.lowerBound();
	int width = upper[0] - lower[0] + 1;
	int height = upper[1] - lower[1] +  1;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			bool shouldBeDilated = process(i, j, size, size);
			if (shouldBeDilated) toReturn.setValue(Point(i,j), 1);
		}
	}
	Image out(domain);
	for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
		out.setValue(*it, toReturn(*it));
	}
	return out;
}

template <typename Image>
Image Morphomaths<Image>::open(int size) {
	myImage = erosion(size);
	Image dilat = dilation(size);
	return dilat;
}

template <typename Image>
Image Morphomaths<Image>::close(int size) {
	myImage = dilation(size);
	Image eros = erosion(size);
	return eros;
}

#endif
