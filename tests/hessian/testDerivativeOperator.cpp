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
 * Functions for testing class DerivativeOperator
 *
 * This file is part of the DGtal library.
 */

#define CATCH_CONFIG_MAIN

#include "catch.hpp"
///////////////////////////////////////////////////////////////////////////////
#include <hessian/DerivativeOperator.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class DerivativeOperator.
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Testing DerivativeOperator") {

    typedef DGtal::DerivativeOperator<double> Derivative;
    typedef Derivative::CoefficientVector Coefficients;

    SECTION("Testing feature computeCoefficients() of DerivativeOperator") {
        Derivative d0(0);
        Derivative d1(1);
        Derivative d2(2);


        Coefficients coeff0 = d0.computeCoefficients();
        Coefficients coeff1 = d1.computeCoefficients();
        Coefficients coeff2 = d2.computeCoefficients();

        Coefficients expectedCoeff1({0.5, 0, -0.5});
        Coefficients expectedCoeff2({1, -2, 1});


        REQUIRE((coeff0.size() == 1));
        REQUIRE((coeff1.size() == 3));
        REQUIRE((coeff2.size() == 3));

        REQUIRE((coeff0[0] == 1));
        REQUIRE((coeff1 == expectedCoeff1));
        REQUIRE((coeff2 == expectedCoeff2));
    }

}

/** @ingroup Tests **/