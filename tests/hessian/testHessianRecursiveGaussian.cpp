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
 * @date 2017/04/05
 *
 * Functions for testing class HessianRecursiveGaussian
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
#include <hessian/itkHessianRecursiveGaussianImageFilter.h>
#include <hessian/HessianRecursiveGaussian.h>
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class HessianRecursiveGaussian.
///////////////////////////////////////////////////////////////////////////////

//TEST_CASE( "Testing 2D HessianRecursiveGaussian" )
//{
//
//    using namespace Z2i;
//    typedef DGtal::ImageContainerBySTLVector<Domain, int> Image;
//    typedef DGtal::ImageContainerByITKImage<Domain, int> DGtalITKImage;
//    typedef DGtal::HessianRecursiveGaussian<Image> Hessian;
//    typedef typename Hessian::OutputImage OutputImage;
//
//
//    typedef itk::Image<int, 2> ITKImage;
//    typedef itk::RecursiveGaussianImageFilter<ITKImage> ITKGaussianOperator;
//    typedef itk::HessianRecursiveGaussianImageFilter<ITKImage> ITKHessian;
//    typedef typename ITKHessian::OutputImageType ITKOutputImage;
//
//
//    Domain domain(Point(0, 0), Point(10, 10));
//    Image image = DGtal::GenericReader<Image>::import(
//            "/home/florent/test_img/Vesselness/DGtal/smallrectangle.pgm");
//    DGtalITKImage image2 = DGtal::GenericReader<DGtalITKImage>::import(
//            "/home/florent/test_img/Vesselness/DGtal/smallrectangle.pgm");
//    ITKImage::Pointer itkImage = image2.getITKImagePointer();
//
////    Image image(domain);
////    for (const Point&  p : domain) {
////        image.setValue(p, 0);
////    }
////    for (int i = 2; i < 8; i++) {
////        for (int j = 2; j < 8; j++) {
////            if (i == 2 && j == 2 ||
////                i == 2 && j == 7 ||
////                i == 7 && j == 2 ||
////                i == 7 && j == 7)
////                continue;
////            image.setValue(Point(i, j), 208);
////        }
////    }
//
//
//    SECTION("Gaussian image") {
//        ITKHessian::Pointer itkHessian = ITKHessian::New();
//        itkHessian->SetInput(itkImage);
//        itkHessian->SetSigma(3);
//        itkHessian->SetNormalizeAcrossScale(true);
//        itkHessian->Update();
//
//        DGtal::trace.info() << "DGtal" << std::endl;
//        ITKOutputImage::Pointer itkOut = itkHessian->GetOutput();
//        Hessian hessian(image, 3.0, false);
//
//        OutputImage hessianImage = hessian.computeHessian();
//        DGtal::trace.info() << hessianImage(Point::zero) << std::endl;
//
//
//
//    }
//}

TEST_CASE("Testing 3D RecursiveGaussianDerivativeOperator") {

    using namespace Z3i;
    typedef DGtal::ImageContainerBySTLVector<Domain, int> Image;
    typedef DGtal::ImageContainerByITKImage<Domain, int> DGtalITKImage;
    typedef DGtal::HessianRecursiveGaussian<Image> Hessian;
    typedef typename Hessian::OutputImage OutputImage;

    typedef itk::Image<int, 3> ITKImage;
    typedef itk::RecursiveGaussianImageFilter<ITKImage> ITKGaussianOperator;
    typedef itk::HessianRecursiveGaussianImageFilter<ITKImage> ITKHessian;
    typedef typename ITKHessian::OutputImageType ITKOutputImage;

    Domain domain(Point(0, 0), Point(10, 10));
    Image image = DGtal::ITKReader<Image>::importITK(
            "/home/florent/test_img/Vesselness/DGtal/smallcylinder.tif");
    DGtalITKImage image2 = DGtal::ITKReader<DGtalITKImage>::importITK(
            "/home/florent/test_img/Vesselness/DGtal/smallcylinder.tif");
    ITKImage::Pointer itkImage = image2.getITKImagePointer();

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
        DGtal::trace.info() << "ITK" << std::endl;

        ITKHessian::Pointer itkHessian = ITKHessian::New();
        itkHessian->SetInput(itkImage);
        itkHessian->SetSigma(2);
        itkHessian->SetNormalizeAcrossScale(true);
        itkHessian->Update();

        DGtal::trace.info() << "DGtal" << std::endl;

        ITKOutputImage::Pointer itkOut = itkHessian->GetOutput();
        Hessian hessian(image, 2.0, true);

        OutputImage hessianImage = hessian.computeHessian();
    }

}

/** @ingroup Tests **/