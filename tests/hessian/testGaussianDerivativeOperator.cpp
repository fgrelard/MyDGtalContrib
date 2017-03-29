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
 * @date 2017/03/28
 *
 * Functions for testing class GaussianDerivativeOperator
 *
 * This file is part of the DGtal library.
 */

#define CATCH_CONFIG_MAIN

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include <catch.hpp>
#include <itkGaussianDerivativeOperator.h>
#include <hessian/GaussianDerivativeOperator.h>
#include "DGtal/helpers/StdDefs.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GaussianDerivativeOperator.
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Testing GaussianDerivativeOperator") {

    typedef std::vector<double> CoefficientVector;
    typedef DGtal::GaussianDerivativeOperator<double> GaussianDer;

    SECTION("Testing feature computeGaussianCoefficients() for order 0 of GaussianDerivativeOperator") {

        GaussianDer gd;
        gd.setOrder(0);
        CoefficientVector coeff0 = gd.computeCoefficients();


        REQUIRE((coeff0.size() == 7));
        REQUIRE((coeff0[3] == Approx(0.466801)));
    }

    SECTION("Testing feature computeCoefficients() for order > 0 of GaussianDerivativeOperator") {

        GaussianDer gd;

        gd.setOrder(1);
        CoefficientVector coeff1 = gd.computeCoefficients();

        gd.setOrder(2);
        CoefficientVector coeff2 = gd.computeCoefficients();

        REQUIRE((coeff1.size() == 7));
        REQUIRE((coeff2.size() == 7));
        REQUIRE((coeff1[3] == 0));
        REQUIRE((coeff2[3] == Approx(-0.516852)));
    }

}

/** @ingroup Tests **/