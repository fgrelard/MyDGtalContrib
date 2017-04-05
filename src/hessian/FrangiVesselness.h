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

#pragma once

/**
 * @file FrangiVesselness.h
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/04/03
 *
 *
 * This file is part of the DGtal library.
 */

#if defined(FrangiVesselness_RECURSES)
#error Recursive header files inclusion detected in FrangiVesselness.h
#else // defined(FrangiVesselness_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FrangiVesselness_RECURSES

#if !defined FrangiVesselness_h
/** Prevents repeated inclusion of headers. */
#define FrangiVesselness_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/math/linalg/EigenDecomposition.h>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class FrangiVesselness
    /**
     * Description of template class 'FrangiVesselness' <p>
     * \brief Aim:
     */
    template <typename THessianImage>
    class FrangiVesselness {
    public:
        typedef THessianImage HessianImage;
        typedef typename HessianImage::Value HessianMatrix;
        typedef typename HessianMatrix::Component HessianValue;
        typedef typename HessianImage::Domain Domain;
        typedef typename Domain::Point Point;
        typedef typename Domain::Space::RealVector RealVector;
        typedef HessianValue VesselnessValue;
        typedef DGtal::ImageContainerBySTLVector<Domain, VesselnessValue> OutputImage;
        typedef DGtal::EigenDecomposition<Domain::dimension, VesselnessValue> LinearAlgebraTool;

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        FrangiVesselness() : myHessianImage(Domain(Point::zero, Point::zero)), myAlpha(0.5), myBeta(1.0),
                             myGamma(10.0) {}

        FrangiVesselness(const HessianImage &aHessianImage) : myHessianImage(aHessianImage), myAlpha(0.5), myBeta(1.0),
                                                              myGamma(10.0) {}

        FrangiVesselness(const HessianImage &aHessianImage, double alpha, double beta, double gamma) : myHessianImage(
                aHessianImage), myAlpha(alpha), myBeta(beta), myGamma(gamma) {}

        /**
         * Destructor.
         */
        ~FrangiVesselness() = default;

        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        FrangiVesselness(const FrangiVesselness &other) = default;

        /**
         * Move constructor.
         * @param other the object to move.
         */
        FrangiVesselness(FrangiVesselness &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        FrangiVesselness &operator=(const FrangiVesselness &other) = default;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        FrangiVesselness &operator=(FrangiVesselness &&other) = delete;

        // ----------------------- Interface --------------------------------------
    public:

        OutputImage computeVesselness();


        /**
         * Writes/Displays the object on an output stream.
         * @param out the output stream where the object is written.
         */
        void selfDisplay(std::ostream &out) const;

        /**
         * Checks the validity/consistency of the object.
         * @return 'true' if the object is valid, 'false' otherwise.
         */
        bool isValid() const;

        HessianImage getImage() const;

        void setImage(const HessianImage &image);

        void setAlpha(double alpha);

        double getAlpha() const;

        void setBeta(double beta);

        double getBeta() const;

        void setGamma(double gamma);

        double getGamma();


        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        HessianImage myHessianImage;
        double myAlpha;
        double myBeta;
        double myGamma;

        // ------------------------- Hidden services ------------------------------
    protected:
        VesselnessValue eigenValueCombination(const RealVector &eigenVal);

        // ------------------------- Internals ------------------------------------
    private:

    }; // end of class FrangiVesselness


    /**
     * Overloads 'operator<<' for displaying objects of class 'FrangiVesselness'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'FrangiVesselness' to write.
     * @return the output stream after the writing.
     */
    template <typename THessianImage>
    std::ostream &
    operator<<(std::ostream &out, const FrangiVesselness<THessianImage> &object) {
        object.selfDisplay(out);
        return out;
    }

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
template <typename THessianImage>
void
DGtal::FrangiVesselness<THessianImage>::
selfDisplay(std::ostream &out) const {
    out << "[FrangiVesselness] image=" << myHessianImage << ", alpha=" << myAlpha << ", beta=" << myBeta << ", gamma="
        << myGamma;
}


/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
template <typename TImage>
bool
DGtal::FrangiVesselness<TImage>::
isValid() const {
    return true;
}

template <typename THessianImage>
typename DGtal::FrangiVesselness<THessianImage>::OutputImage
DGtal::FrangiVesselness<THessianImage>::
computeVesselness() {
    OutputImage frangiImage(myHessianImage.domain());
    for (const Point &p : myHessianImage.domain()) {
        HessianMatrix hessian = myHessianImage(p);
        HessianMatrix eigvec;
        RealVector eigval;
        LinearAlgebraTool::getEigenDecomposition(hessian, eigvec, eigval);
        VesselnessValue frangiValue = eigenValueCombination(eigval);
        frangiImage.setValue(p, frangiValue);
    }
    return frangiImage;
}

template <typename THessianImage>
typename DGtal::FrangiVesselness<THessianImage>::VesselnessValue
DGtal::FrangiVesselness<THessianImage>::
eigenValueCombination(const RealVector &eigval) {
    RealVector sortedEigVal = eigval;
    std::sort(sortedEigVal.begin(), sortedEigVal.end(), [&](HessianValue h1, HessianValue h2) {
        return std::abs(h1) < std::abs(h2);
    });
    HessianValue l1 = std::abs(sortedEigVal[0]);
    HessianValue l2 = std::abs(sortedEigVal[1]);

    switch (Domain::dimension) {
        case 2: {
            if (sortedEigVal[1] > DGtal::NumberTraits<HessianValue>::ZERO)
                return DGtal::NumberTraits<VesselnessValue>::ZERO;
            VesselnessValue rb = l1 / l2;
            VesselnessValue s = std::sqrt(l1 * l1 + l2 * l2);
            VesselnessValue frangi = std::exp(-(rb * rb) / (2.0 * myBeta * myBeta));
            frangi *= 1.0 - std::exp(-(s * s) / (2.0 * myGamma * myGamma));
            return frangi;
        }

        case 3: {
            if (sortedEigVal[1] > DGtal::NumberTraits<HessianValue>::ZERO ||
                sortedEigVal[2] > DGtal::NumberTraits<HessianValue>::ZERO)
                return DGtal::NumberTraits<VesselnessValue>::ZERO;
            HessianValue l3 = std::abs(sortedEigVal[2]);
            VesselnessValue ra = l2 / l3;
            VesselnessValue rb = l1 / std::sqrt(l2 * l3);
            VesselnessValue s = std::sqrt(l1 * l1 + l2 * l2 + l3 * l3);
            VesselnessValue frangi = 1.0 - std::exp(-(ra * ra) / (2.0 * myAlpha * myAlpha));
            frangi *= std::exp(-(rb * rb) / (2.0 * myBeta * myBeta));
            frangi *= 1.0 - std::exp(-(s * s) / (2.0 * myGamma * myGamma));
            return frangi;
        }
        default:
            return DGtal::NumberTraits<VesselnessValue>::ZERO;
    }

}

template <typename THessianImage>
void
DGtal::FrangiVesselness<THessianImage>::
setImage(const HessianImage &image) {
    myHessianImage = image;
}

template <typename THessianImage>
typename DGtal::FrangiVesselness<THessianImage>::HessianImage
DGtal::FrangiVesselness<THessianImage>::
getImage() const {
    return myHessianImage;
}

template <typename THessianImage>
void
DGtal::FrangiVesselness<THessianImage>::
setAlpha(double alpha) {
    myAlpha = alpha;
}

template <typename THessianImage>
double
DGtal::FrangiVesselness<THessianImage>::
getAlpha() const {
    return myAlpha;
}

template <typename THessianImage>
void
DGtal::FrangiVesselness<THessianImage>::
setBeta(double beta) {
    myBeta = beta;
}

template <typename THessianImage>
double
DGtal::FrangiVesselness<THessianImage>::
getBeta() const {
    return myBeta;
}

template <typename THessianImage>
void
DGtal::FrangiVesselness<THessianImage>::
setGamma(double gamma) {
    myGamma = gamma;
}

template <typename THessianImage>
double
DGtal::FrangiVesselness<THessianImage>::
getGamma() {
    return myGamma;
}

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FrangiVesselness_h

#undef FrangiVesselness_RECURSES
#endif // else defined(FrangiVesselness_RECURSES)