/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/


/**
 * @file
 * @ingroup Tests
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/04/04
 *
 * Functions for testing class RecursiveGaussianDerivativeOperator
 *
 * This file is part of the DGtal library.
 */

#define CATCH_CONFIG_MAIN
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include <catch.hpp>
#include <DGtal/io/readers/GenericReader.h>
#include <hessian/RecursiveGaussianDerivativeOperator.h>
#include <itkRecursiveGaussianImageFilter.h>
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class RecursiveGaussianDerivativeOperator.
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Testing 2D RecursiveGaussianDerivativeOperator") {

    using namespace Z2i;
    typedef DGtal::ImageContainerBySTLVector<Domain, int> Image;
    typedef DGtal::ImageContainerByITKImage<Domain, int> DGtalITKImage;
    typedef DGtal::RecursiveGaussianDerivativeOperator<DGtalITKImage> GaussianOperator;
    typedef typename GaussianOperator::OutputImage OutputImage;


    typedef itk::Image<int, 2> ITKImage;
    typedef itk::RecursiveGaussianImageFilter<ITKImage> ITKGaussianOperator;
    typedef typename ITKGaussianOperator::OutputImageType ITKOutputImage;

    Domain domain(Point(0, 0), Point(10, 10));
    DGtalITKImage image = DGtal::GenericReader<DGtalITKImage>::import(
            "/home/florent/test_img/Vesselness/DGtal/smallcircle.pgm");
    ITKImage::Pointer itkImage = image.getITKImagePointer();

//    Image image(domain);
//    for (const Point&  p : domain) {
//        image.setValue(p, 0);
//    }
//    for (int i = 2; i < 8; i++) {
//        for (int j = 2; j < 8; j++) {
//            if (i == 2 && j == 2 ||
//                i == 2 && j == 7 ||
//                i == 7 && j == 2 ||
//                i == 7 && j == 7)
//                continue;
//            image.setValue(Point(i, j), 208);
//        }
//    }


    SECTION("Gaussian image") {
        ITKGaussianOperator::Pointer itkGaussianOperator = ITKGaussianOperator::New();
        itkGaussianOperator->SetInput(itkImage);
        itkGaussianOperator->SetFirstOrder();
        itkGaussianOperator->SetSigma(3.0);
        itkGaussianOperator->SetNormalizeAcrossScale(true);
        itkGaussianOperator->Update();
        GaussianOperator gaussianOperator(image, 3.0, 1);
        OutputImage out = gaussianOperator.gaussianImage();
        ITKOutputImage::Pointer itkOut = itkGaussianOperator->GetOutput();
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                typename ITKOutputImage::IndexType index;
                index[0] = i;
                index[1] = j;

                DGtal::trace.info() << out(Point(i, j)) << " ";
            }
            DGtal::trace.info() << std::endl;
        }

        DGtal::trace.info() << "ITK" << std::endl;
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                typename ITKOutputImage::IndexType index;
                index[0] = i;
                index[1] = j;

                DGtal::trace.info() << itkOut->GetPixel(index) << " ";
            }
            DGtal::trace.info() << std::endl;
        }

    }

}

/** @ingroup Tests **/