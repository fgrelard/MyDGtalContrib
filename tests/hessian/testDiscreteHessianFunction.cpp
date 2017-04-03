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
 * @date 2017/03/31
 *
 * Functions for testing class DiscreteHessianFunction
 *
 * This file is part of the DGtal library.
 */


#define CATCH_CONFIG_MAIN

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include <catch.hpp>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <hessian/DiscreteHessianFunction.h>
#include "DGtal/helpers/StdDefs.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class DiscreteHessianFunction.
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Testing DiscreteHessianFunction 2D") {
    typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D;
    typedef DGtal::DiscreteHessianFunction<Image2D> Hessian2D;

    Z2i::Domain domain(Z2i::Point::zero, Z2i::Point(10, 10));
    Image2D image2D(domain);
    for (const Z2i::Point &p : domain) {
        if (p[0] == 5)
            image2D.setValue(p, 120);
    }
    Hessian2D hessian2D(image2D, 1.0);

    SECTION("Initial hessian") {
        typename Hessian2D::OutputImage im = hessian2D.computeHessian();
        typename Hessian2D::HessianMatrix mat = *im.begin();
        Z2i::Point pointTube(5, 2);
        typename Hessian2D::HessianMatrix matTube = im(pointTube);

        REQUIRE((mat(0, 0) == 0.0));
        REQUIRE((matTube(0, 0) == Approx(-61.5153) && matTube(0, 1) == Approx(0) && matTube(1, 1) == Approx(-2.34578)));
    }

    SECTION("Changing sigma") {
        Z2i::Point pointTube(5, 2);
        typename Hessian2D::OutputImage imBefore = hessian2D.computeHessian();
        typename Hessian2D::HessianMatrix matBefore = imBefore(pointTube);

        hessian2D.setSigma(0.5);
        typename Hessian2D::OutputImage imAfter = hessian2D.computeHessian();
        typename Hessian2D::HessianMatrix matAfter = imAfter(pointTube);
        REQUIRE((matBefore(0, 0) != matAfter(0, 0) || matBefore(0, 1) != matAfter(0, 1) ||
                 matBefore(1, 0) != matAfter(1, 1)));

    }

    SECTION("Evaluating hessian at point inside domain") {
        Z2i::Point pointTube(5, 2);
        typename Hessian2D::HessianMatrix matInside = hessian2D.hessianAtPoint(pointTube);
        REQUIRE((matInside(0, 0) == Approx(-61.5153) && matInside(0, 1) == Approx(0) &&
                 matInside(1, 1) == Approx(-2.34578)));
    }

    SECTION("Evaluating hessian at point outside domain") {
        Z2i::Point pointTube(11, 2);
        typename Hessian2D::HessianMatrix matOutside = hessian2D.hessianAtPoint(pointTube);
        REQUIRE((matOutside(0, 0) == Approx(0) && matOutside(0, 1) == Approx(0) && matOutside(1, 1) == Approx(0)));
    }


}

TEST_CASE("Testing DiscreteHessianFunction 3D") {

    typedef DGtal::ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3D;
    typedef DGtal::DiscreteHessianFunction<Image3D>
            Hessian3D;


    Z3i::Domain domain(Z3i::Point::zero, Z3i::Point::diagonal(10));
    Image3D image3D(domain);
    for (const Z3i::Point &p : domain) {
        image3D.setValue(p, 0);
    }
    for (int i = 2; i < 8; i++) {
        for (int j = 2; j < 8; j++) {
            for (int k = 2; k < 8; k++) {
                if (i == 2 && j == 2 ||
                    i == 2 && j == 7 ||
                    i == 7 && j == 2 ||
                    i == 7 && j == 7)
                    continue;
                image3D.setValue(Z3i::Point(i, j, k), 208);
            }
        }
    }


    Hessian3D hessian3D(image3D, 1.0);

    SECTION("Initial hessian") {
        typename Hessian3D::OutputImage im = hessian3D.computeHessian();
        Z3i::Point pointTube(6, 6, 2);
        typename Hessian3D::HessianMatrix matTube = im(pointTube);
        REQUIRE((matTube(0, 0) == Approx(-25.9277) &&
                 matTube(0, 1) == Approx(-4.38875) &&
                 matTube(0, 2) == Approx(-10.7856) &&
                 matTube(1, 1) == Approx(-25.9277) &&
                 matTube(1, 2) == Approx(-10.7856) &&
                 matTube(2, 2) == Approx(-45.3415)));
    }

    SECTION("Changing sigma") {
        Z3i::Point pointTube(5, 2);
        typename Hessian3D::OutputImage imBefore = hessian3D.computeHessian();
        typename Hessian3D::HessianMatrix matBefore = imBefore(pointTube);

        hessian3D.setSigma(0.5);
        typename Hessian3D::OutputImage imAfter = hessian3D.computeHessian();
        typename Hessian3D::HessianMatrix matAfter = imAfter(pointTube);
        REQUIRE((matBefore(0, 0) != matAfter(0, 0) || matBefore(0, 1) != matAfter(0, 1) ||
                 matBefore(1, 0) != matAfter(1, 1)));

    }

    SECTION("Evaluating hessian at point inside domain") {

        Z3i::Point pointTube(6, 6, 2);
        typename Hessian3D::HessianMatrix matTube = hessian3D.hessianAtPoint(pointTube);
        REQUIRE((matTube(0, 0) == Approx(-25.9277) &&
                 matTube(0, 1) == Approx(-4.38875) &&
                 matTube(0, 2) == Approx(-10.7856) &&
                 matTube(1, 1) == Approx(-25.9277) &&
                 matTube(1, 2) == Approx(-10.7856) &&
                 matTube(2, 2) == Approx(-45.3415)));


    }

    SECTION("Evaluating hessian at point outside domain") {
        Z3i::Point pointTube(11, 2, 2);
        typename Hessian3D::HessianMatrix matOutside = hessian3D.hessianAtPoint(pointTube);
        REQUIRE((matOutside(0, 0) == Approx(0) && matOutside(0, 1) == Approx(0) && matOutside(0, 2) == Approx(0)));
    }

}

/** @ingroup Tests **/