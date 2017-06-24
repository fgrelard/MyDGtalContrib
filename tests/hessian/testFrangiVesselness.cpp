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
 * Functions for testing class FrangiVesselness
 *
 * This file is part of the DGtal library.
 */

#define CATCH_CONFIG_MAIN
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include <catch.hpp>
#include "DGtal/helpers/StdDefs.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class FrangiVesselness.
///////////////////////////////////////////////////////////////////////////////

TEST_CASE( "Testing 2D FrangiVesselness" )
{

    using namespace Z2i;
    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image;
    typedef DGtal::DiscreteHessianFunction<Image> Hessian;
    typedef typename Hessian::OutputImage HessianImage;
    typedef DGtal::FrangiVesselness<HessianImage> Vesselness;

    SECTION("Testing feature ZZZ of FrangiVesselness") {
        REQUIRE( (a == b) );
    }

}

/** @ingroup Tests **/